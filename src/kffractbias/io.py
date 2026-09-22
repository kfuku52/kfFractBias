from __future__ import annotations

import gzip
import hashlib
import re
from collections import defaultdict
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import TextIO
from urllib.parse import unquote


@dataclass(frozen=True)
class Gene:
    seqid: str
    start: int
    end: int
    gene_id: str
    strand: str = "."


@dataclass(frozen=True)
class AnnotationMapping:
    feature: str
    attribute: str
    genes: tuple[Gene, ...]
    fasta_gene_count: int
    matched_gene_count: int
    locus_by_id: dict[str, str] = field(default_factory=dict)
    isoform_policy: str = "error"
    collapsed_isoform_count: int = 0
    unresolved_locus_count: int = 0

    def metadata(self) -> dict[str, str | int]:
        return {
            "feature": self.feature,
            "attribute": self.attribute,
            "fasta_gene_count": self.fasta_gene_count,
            "matched_gene_count": self.matched_gene_count,
            "selected_gene_count": len(self.genes),
            "isoform_policy": self.isoform_policy,
            "collapsed_isoform_count": self.collapsed_isoform_count,
            "unresolved_locus_count": self.unresolved_locus_count,
            "counting_unit": "mapped_identifier"
            if self.isoform_policy == "all"
            else (
                "representative_per_known_locus"
                if self.unresolved_locus_count
                else "representative_per_locus"
            ),
        }


@dataclass(frozen=True)
class SyntenyParseResult:
    pairs: tuple[tuple[str, str], ...]
    record_count: int
    duplicate_pair_count: int


@dataclass(frozen=True)
class _GffRecord:
    seqid: str
    feature: str
    start: int
    end: int
    strand: str
    attributes: dict[str, tuple[str, ...]]
    line_number: int


FEATURE_PRIORITY = ("mRNA", "transcript", "gene", "CDS")
ATTRIBUTE_PRIORITY = (
    "ID",
    "Name",
    "transcript_id",
    "gene_id",
    "protein_id",
    "locus_tag",
    "Parent",
    "coge_fid",
)


@contextmanager
def open_text(path: str | Path) -> Iterator[TextIO]:
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            yield handle
    else:
        with path.open("r", encoding="utf-8") as handle:
            yield handle


def natural_key(value: str) -> tuple[object, ...]:
    parts = tuple(
        int(token) if token.isdecimal() else token.lower() for token in re.split(r"(\d+)", value)
    )
    return parts, value


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_gene_identifier(value: str, format_name: str, location: str) -> None:
    if not value:
        raise ValueError(f"Empty {format_name} gene identifier at {location}")
    if (
        value.startswith("#")
        or ";" in value
        or "||" in value
        or any(c.isspace() or ord(c) < 32 or ord(c) == 127 for c in value)
    ):
        raise ValueError(
            f"Invalid {format_name} gene identifier {value!r} at {location}; "
            "identifiers must not start with # or contain whitespace, control characters, ; or ||"
        )


def validate_disjoint_identifiers(
    target_ids: Iterable[str], query_ids: Iterable[str], label: str
) -> None:
    overlapping = set(target_ids).intersection(query_ids)
    if overlapping:
        examples = ", ".join(sorted(overlapping, key=natural_key)[:10])
        raise ValueError(
            f"Pairwise target and query {label} identifiers must be disjoint; "
            f"found {len(overlapping)} overlapping identifier(s), including: {examples}"
        )


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str, str]]:
    """Validate nucleotide FASTA while retaining at most one sequence record."""
    identifiers: set[str] = set()
    identifier = header = ""
    header_line = 0
    sequence: list[str] = []
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            if line.startswith(">"):
                if identifier:
                    if not sequence:
                        raise ValueError(
                            f"Empty FASTA sequence for {identifier!r} at {path}:{header_line}"
                        )
                    yield identifier, header, "".join(sequence)
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA identifier at {path}:{line_number}")
                identifier = header.split(None, 1)[0]
                validate_gene_identifier(identifier, "FASTA", f"{path}:{line_number}")
                if identifier in identifiers:
                    raise ValueError(
                        f"Duplicate FASTA identifier {identifier!r} at {path}:{line_number}"
                    )
                identifiers.add(identifier)
                header_line = line_number
                sequence = []
            else:
                if not identifier:
                    raise ValueError(f"Sequence before first FASTA header at {path}:{line_number}")
                bases = "".join(line.split())
                if re.search(r"[^ACGTURYSWKMBDHVNacgturyswkmbdhvn]", bases):
                    raise ValueError(f"Invalid nucleotide FASTA sequence at {path}:{line_number}")
                sequence.append(bases)
    if not identifiers:
        raise ValueError(f"No FASTA records found in {path}")
    if not sequence:
        raise ValueError(f"Empty FASTA sequence for {identifier!r} at {path}:{header_line}")
    yield identifier, header, "".join(sequence)


def read_fasta_ids(path: str | Path) -> set[str]:
    return {identifier for identifier, _header, _sequence in iter_fasta(path)}


def select_isoforms(
    mapping: AnnotationMapping, lengths: dict[str, int], policy: str
) -> AnnotationMapping:
    if policy not in {"error", "longest", "all"}:
        raise ValueError("isoform policy must be 'error', 'longest', or 'all'")
    loci: dict[str, list[Gene]] = defaultdict(list)
    for gene in mapping.genes:
        loci[mapping.locus_by_id.get(gene.gene_id, gene.gene_id)].append(gene)
    duplicated = {locus: genes for locus, genes in loci.items() if len(genes) > 1}
    if duplicated and policy == "error":
        examples = ", ".join(sorted(duplicated, key=natural_key)[:5])
        raise ValueError(
            f"Multiple CDS isoforms map to {len(duplicated)} gene loci ({examples}); use --isoform-policy longest to select one CDS per locus, or all to count identifiers separately"
        )
    selected = mapping.genes
    if policy == "longest":
        for locus, genes in duplicated.items():
            if (
                len({gene.seqid for gene in genes}) > 1
                or len({gene.strand for gene in genes} - {"."}) > 1
            ):
                raise ValueError(
                    f"Isoforms of locus {locus!r} occur on incompatible sequences or strands"
                )
        ids = {
            min(genes, key=lambda gene: (-lengths[gene.gene_id], natural_key(gene.gene_id))).gene_id
            for genes in loci.values()
        }
        selected = tuple(gene for gene in selected if gene.gene_id in ids)
    return replace(
        mapping,
        genes=selected,
        isoform_policy=policy,
        collapsed_isoform_count=len(mapping.genes) - len(selected),
    )


def parse_attributes(value: str) -> dict[str, tuple[str, ...]]:
    parsed: dict[str, list[str]] = defaultdict(list)
    for attribute_field in value.strip().strip(";").split(";"):
        attribute_field = attribute_field.strip()
        if not attribute_field:
            continue
        if "=" in attribute_field:
            key, raw = attribute_field.split("=", 1)
        else:
            match = re.match(r"([^\s]+)\s+[\"']?(.*?)[\"']?$", attribute_field)
            if not match:
                continue
            key, raw = match.groups()
        for item in raw.split(","):
            item = unquote(item.strip().strip("\"'"))
            if item:
                parsed[key].append(item)
    return {key: tuple(values) for key, values in parsed.items()}


def _validate_seqid(value: str, format_name: str, location: str) -> None:
    if not value.strip():
        raise ValueError(f"Empty {format_name} sequence identifier at {location}")
    if value.startswith("#") or any(
        character.isspace() or ord(character) < 32 or ord(character) == 127 for character in value
    ):
        raise ValueError(f"Invalid {format_name} sequence identifier at {location}")


def _parse_strand(value: str, format_name: str, location: str) -> str:
    if value not in {"+", "-", ".", "?"}:
        raise ValueError(f"Invalid {format_name} strand {value!r} at {location}")
    return "." if value == "?" else value


def _iter_gff(path: str | Path) -> Iterator[_GffRecord]:
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9:
                raise ValueError(
                    f"Expected 9 GFF columns at {path}:{line_number}; found {len(columns)}"
                )
            _validate_seqid(columns[0], "GFF", f"{path}:{line_number}")
            strand = _parse_strand(columns[6], "GFF", f"{path}:{line_number}")
            try:
                start = int(columns[3]) - 1
                end = int(columns[4])
            except ValueError as exc:
                raise ValueError(f"Invalid GFF coordinates at {path}:{line_number}") from exc
            if start < 0 or end <= start:
                raise ValueError(f"Invalid GFF interval at {path}:{line_number}")
            yield _GffRecord(
                columns[0],
                columns[2],
                start,
                end,
                strand,
                parse_attributes(columns[8]),
                line_number,
            )


def _select_gff_mapping(
    gff_path: str | Path,
    fasta_ids: set[str],
    features: tuple[str, ...],
    attributes: tuple[str, ...],
) -> tuple[str, str, set[str]]:
    matches: dict[tuple[str, str], set[str]] = {
        (candidate_feature, candidate_attribute): set()
        for candidate_feature in features
        for candidate_attribute in attributes
    }
    observed_features: set[str] = set()
    observed_attributes: set[str] = set()

    for row in _iter_gff(gff_path):
        observed_features.add(row.feature)
        observed_attributes.update(row.attributes)
        if row.feature not in features:
            continue
        for candidate_attribute in attributes:
            for candidate_id in row.attributes.get(candidate_attribute, ()):
                if candidate_id in fasta_ids:
                    matches[(row.feature, candidate_attribute)].add(candidate_id)

    selected_pair, selected_matches = max(
        matches.items(),
        key=lambda item: (
            len(item[1]),
            -features.index(item[0][0]),
            -attributes.index(item[0][1]),
        ),
    )
    if not selected_matches:
        observed_feature_text = ", ".join(sorted(observed_features, key=natural_key)) or "none"
        observed_attribute_text = ", ".join(sorted(observed_attributes, key=natural_key)) or "none"
        raise ValueError(
            "No GFF identifiers matched the CDS FASTA. "
            f"Observed features: {observed_feature_text}. Observed attributes: {observed_attribute_text}."
        )

    return selected_pair[0], selected_pair[1], selected_matches


def _mapped_gff_intervals(
    gff_path: str | Path,
    selected_feature: str,
    selected_attribute: str,
    selected_matches: set[str],
) -> tuple[Gene, ...]:
    intervals: dict[str, Gene] = {}
    for row in _iter_gff(gff_path):
        if row.feature != selected_feature:
            continue
        for gene_id in row.attributes.get(selected_attribute, ()):
            if gene_id not in selected_matches:
                continue
            validate_gene_identifier(gene_id, "GFF", f"{gff_path}:{row.line_number}")
            previous = intervals.get(gene_id)
            if previous is None:
                intervals[gene_id] = Gene(row.seqid, row.start, row.end, gene_id, row.strand)
                continue
            location = f"{gff_path}:{row.line_number}"
            if previous.seqid != row.seqid:
                raise ValueError(
                    f"GFF identifier {gene_id!r} occurs on multiple sequences at {location}"
                )
            known_strands = {previous.strand, row.strand} - {"."}
            if len(known_strands) > 1:
                raise ValueError(
                    f"GFF identifier {gene_id!r} has conflicting strands at {location}"
                )
            intervals[gene_id] = Gene(
                row.seqid,
                min(previous.start, row.start),
                max(previous.end, row.end),
                gene_id,
                next(iter(known_strands), "."),
            )

    return tuple(
        sorted(
            intervals.values(),
            key=lambda gene: (natural_key(gene.seqid), gene.start, gene.end, gene.gene_id),
        )
    )


def _annotation_loci(
    gff_path: str | Path, feature: str, attribute: str, selected: set[str]
) -> tuple[dict[str, str], int]:
    parents: dict[str, set[str]] = defaultdict(set)
    candidates: dict[str, set[str]] = defaultdict(set)
    known_loci: set[str] = set()
    for row in _iter_gff(gff_path):
        row_feature, attrs = row.feature, row.attributes
        row_ids = set(attrs.get("ID", ())) | set(attrs.get("transcript_id", ()))
        ancestors = set(attrs.get("gene_id", ())) or set(attrs.get("Parent", ()))
        known_loci.update(attrs.get("gene_id", ()))
        if row_feature == "gene":
            known_loci.update(row_ids)
        elif row_feature in {"mRNA", "transcript"}:
            known_loci.update(ancestors)
        if row_feature != "gene":
            for identifier in row_ids:
                parents[identifier].update(ancestors - {identifier})
        if row_feature != feature:
            continue
        for identifier in set(attrs.get(attribute, ())) & selected:
            if row_feature == "gene":
                candidates[identifier].update(
                    attrs.get("ID", ()) or attrs.get("gene_id", ()) or (identifier,)
                )
                known_loci.update(candidates[identifier])
            elif attribute == "Parent" and not attrs.get("gene_id"):
                # A shared CDS segment can name several transcripts. Resolve
                # each mapped transcript separately, not all of its siblings.
                candidates[identifier].add(identifier)
            else:
                candidates[identifier].update(ancestors)

    def roots(identifier: str, seen: frozenset[str] = frozenset()) -> set[str]:
        if identifier in seen:
            raise ValueError(f"Cyclic GFF parent relationship for {identifier!r} in {gff_path}")
        if not parents.get(identifier):
            return {identifier}
        return set().union(*(roots(parent, seen | {identifier}) for parent in parents[identifier]))

    loci: dict[str, str] = {}
    unresolved = 0
    for identifier in sorted(selected, key=natural_key):
        ancestors = candidates.get(identifier, set())
        if not ancestors:
            unresolved += 1
            loci[identifier] = identifier
            continue
        resolved = set().union(*(roots(parent) for parent in ancestors))
        if len(resolved) != 1:
            raise ValueError(
                f"GFF identifier {identifier!r} maps to multiple gene loci: {', '.join(sorted(resolved))}"
            )
        locus = next(iter(resolved))
        loci[identifier] = locus
        if locus not in known_loci:
            unresolved += 1
    return loci, unresolved


def annotation_to_genes(
    gff_path: str | Path,
    fasta_ids: set[str],
    *,
    feature: str | None = None,
    attribute: str | None = None,
) -> AnnotationMapping:
    if (feature is None) != (attribute is None):
        raise ValueError("GFF feature and attribute must be specified together")

    features = (feature,) if feature else FEATURE_PRIORITY
    attributes = (attribute,) if attribute else ATTRIBUTE_PRIORITY
    selected_feature, selected_attribute, selected_matches = _select_gff_mapping(
        gff_path,
        fasta_ids,
        features,
        attributes,
    )
    genes = _mapped_gff_intervals(
        gff_path,
        selected_feature,
        selected_attribute,
        selected_matches,
    )
    loci, unresolved = _annotation_loci(
        gff_path, selected_feature, selected_attribute, selected_matches
    )
    return AnnotationMapping(
        feature=selected_feature,
        attribute=selected_attribute,
        genes=genes,
        fasta_gene_count=len(fasta_ids),
        matched_gene_count=len(genes),
        locus_by_id=loci,
        unresolved_locus_count=unresolved,
    )


def write_bed(genes: Iterable[Gene], path: str | Path) -> None:
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        handle.writelines(
            f"{gene.seqid}\t{gene.start}\t{gene.end}\t{gene.gene_id}\t0\t{gene.strand}\n"
            for gene in genes
        )


def read_bed(path: str | Path) -> tuple[Gene, ...]:
    genes: list[Gene] = []
    seen: set[str] = set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 4:
                raise ValueError(f"Expected at least 4 BED columns at {path}:{line_number}")
            gene_id = columns[3]
            _validate_seqid(columns[0], "BED", f"{path}:{line_number}")
            validate_gene_identifier(gene_id, "BED", f"{path}:{line_number}")
            if gene_id in seen:
                raise ValueError(f"Duplicate BED identifier {gene_id!r} in {path}")
            seen.add(gene_id)
            try:
                start, end = int(columns[1]), int(columns[2])
            except ValueError as exc:
                raise ValueError(f"Invalid BED coordinates at {path}:{line_number}") from exc
            if start < 0 or end <= start:
                raise ValueError(f"Invalid BED interval at {path}:{line_number}")
            strand = _parse_strand(
                columns[5] if len(columns) >= 6 else ".", "BED", f"{path}:{line_number}"
            )
            genes.append(Gene(columns[0], start, end, gene_id, strand))
    if not genes:
        raise ValueError(f"No BED records found in {path}")
    return tuple(
        sorted(
            genes, key=lambda gene: (natural_key(gene.seqid), gene.start, gene.end, gene.gene_id)
        )
    )


def _unique_member(candidates: Iterable[str], identifiers: set[str]) -> str | None:
    matches = tuple(
        dict.fromkeys(candidate for candidate in candidates if candidate in identifiers)
    )
    if len(matches) > 1:
        raise ValueError(f"Multiple candidate identifiers match the same BED: {', '.join(matches)}")
    return matches[0] if matches else None


def _parse_jcvi_pairs(
    path: str | Path,
    target_ids: set[str],
    query_ids: set[str],
    *,
    allow_ambiguous_orientation: bool,
) -> SyntenyParseResult:
    pairs: set[tuple[str, str]] = set()
    record_count = 0
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            record_count += 1
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 2:
                raise ValueError(f"Expected at least 2 JCVI columns at {path}:{line_number}")
            left, right = columns[0], columns[1]
            forward = left in target_ids and right in query_ids
            reverse = right in target_ids and left in query_ids
            if forward and reverse and not allow_ambiguous_orientation:
                raise ValueError(
                    f"Ambiguous JCVI pair orientation at {path}:{line_number}; "
                    "target and query BED identifiers must not overlap"
                )
            if forward:
                pairs.add((left, right))
            elif reverse:
                pairs.add((right, left))
            else:
                raise ValueError(
                    f"JCVI pair at {path}:{line_number} does not match one target and one query BED identifier: "
                    f"{left!r}, {right!r}"
                )
    if not pairs:
        raise ValueError(f"No JCVI anchor pairs found in {path}")
    ordered_pairs = tuple(sorted(pairs))
    return SyntenyParseResult(
        pairs=ordered_pairs,
        record_count=record_count,
        duplicate_pair_count=record_count - len(ordered_pairs),
    )


def read_jcvi_pairs(
    path: str | Path, target_ids: set[str], query_ids: set[str]
) -> tuple[tuple[str, str], ...]:
    return _parse_jcvi_pairs(
        path,
        target_ids,
        query_ids,
        allow_ambiguous_orientation=False,
    ).pairs


def _parse_synmap_pairs(
    path: str | Path,
    target_ids: set[str],
    query_ids: set[str],
    *,
    allow_ambiguous_orientation: bool,
) -> SyntenyParseResult:
    pairs: set[tuple[str, str]] = set()
    record_count = 0
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            record_count += 1
            columns = line.replace("||", "\t").rstrip("\n").split("\t")
            if len(columns) < 20:
                raise ValueError(
                    f"Expected at least 20 expanded SynMap columns at {path}:{line_number}"
                )
            left_candidates = (columns[4], columns[7])
            right_candidates = (columns[16], columns[19])
            try:
                left_target = _unique_member(left_candidates, target_ids)
                right_query = _unique_member(right_candidates, query_ids)
                right_target = _unique_member(right_candidates, target_ids)
                left_query = _unique_member(left_candidates, query_ids)
            except ValueError as exc:
                raise ValueError(
                    f"Ambiguous SynMap identifiers at {path}:{line_number}: {exc}"
                ) from exc
            if (
                left_target is not None
                and right_query is not None
                and right_target is not None
                and left_query is not None
                and not allow_ambiguous_orientation
            ):
                raise ValueError(
                    f"Ambiguous SynMap pair orientation at {path}:{line_number}; "
                    "target and query BED identifiers must not overlap"
                )
            if left_target is not None and right_query is not None:
                pairs.add((left_target, right_query))
            elif right_target is not None and left_query is not None:
                pairs.add((right_target, left_query))
            else:
                raise ValueError(
                    f"SynMap pair at {path}:{line_number} does not match one target and one query BED identifier"
                )
    if not pairs:
        raise ValueError(f"No SynMap pairs found in {path}")
    ordered_pairs = tuple(sorted(pairs))
    return SyntenyParseResult(
        pairs=ordered_pairs,
        record_count=record_count,
        duplicate_pair_count=record_count - len(ordered_pairs),
    )


def read_synmap_pairs(
    path: str | Path, target_ids: set[str], query_ids: set[str]
) -> tuple[tuple[str, str], ...]:
    return _parse_synmap_pairs(
        path,
        target_ids,
        query_ids,
        allow_ambiguous_orientation=False,
    ).pairs


def detect_synteny_format(path: str | Path) -> str:
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            return "synmap" if "||" in line else "jcvi"
    raise ValueError(f"No synteny records found in {path}")


def read_synteny_pairs(
    path: str | Path,
    format_name: str,
    target_ids: set[str],
    query_ids: set[str],
) -> tuple[tuple[str, str], ...]:
    return parse_synteny_pairs(path, format_name, target_ids, query_ids).pairs


def parse_synteny_pairs(
    path: str | Path,
    format_name: str,
    target_ids: set[str],
    query_ids: set[str],
    *,
    allow_ambiguous_orientation: bool = False,
) -> SyntenyParseResult:
    if format_name == "auto":
        format_name = detect_synteny_format(path)
    if format_name == "jcvi":
        return _parse_jcvi_pairs(
            path,
            target_ids,
            query_ids,
            allow_ambiguous_orientation=allow_ambiguous_orientation,
        )
    if format_name == "synmap":
        return _parse_synmap_pairs(
            path,
            target_ids,
            query_ids,
            allow_ambiguous_orientation=allow_ambiguous_orientation,
        )
    raise ValueError(f"Unsupported synteny format: {format_name}")

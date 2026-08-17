from __future__ import annotations

import gzip
import hashlib
import re
from collections import defaultdict
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
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
    return tuple(int(token) if token.isdigit() else token.lower() for token in re.split(r"(\d+)", value))


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fasta_ids(path: str | Path) -> set[str]:
    identifiers: set[str] = set()
    with open_text(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            identifier = line[1:].strip().split(None, 1)[0]
            if not identifier:
                raise ValueError(f"Empty FASTA identifier in {path}")
            if identifier in identifiers:
                raise ValueError(f"Duplicate FASTA identifier {identifier!r} in {path}")
            identifiers.add(identifier)
    if not identifiers:
        raise ValueError(f"No FASTA records found in {path}")
    return identifiers


def parse_attributes(value: str) -> dict[str, tuple[str, ...]]:
    parsed: dict[str, list[str]] = defaultdict(list)
    for field in value.strip().strip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, raw = field.split("=", 1)
        else:
            match = re.match(r"([^\s]+)\s+[\"']?(.*?)[\"']?$", field)
            if not match:
                continue
            key, raw = match.groups()
        for item in raw.split(","):
            item = unquote(item.strip().strip('"\''))
            if item:
                parsed[key].append(item)
    return {key: tuple(values) for key, values in parsed.items()}


def _iter_gff(path: str | Path) -> Iterator[tuple[str, str, int, int, str, dict[str, tuple[str, ...]]]]:
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9:
                raise ValueError(f"Expected 9 GFF columns at {path}:{line_number}; found {len(columns)}")
            try:
                start = int(columns[3]) - 1
                end = int(columns[4])
            except ValueError as exc:
                raise ValueError(f"Invalid GFF coordinates at {path}:{line_number}") from exc
            if start < 0 or end <= start:
                raise ValueError(f"Invalid GFF interval at {path}:{line_number}")
            yield columns[0], columns[2], start, end, columns[6], parse_attributes(columns[8])


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
    matches: dict[tuple[str, str], set[str]] = {
        (candidate_feature, candidate_attribute): set()
        for candidate_feature in features
        for candidate_attribute in attributes
    }
    observed_features: set[str] = set()
    observed_attributes: set[str] = set()

    for _seqid, row_feature, _start, _end, _strand, attrs in _iter_gff(gff_path):
        observed_features.add(row_feature)
        observed_attributes.update(attrs)
        if row_feature not in features:
            continue
        for candidate_attribute in attributes:
            for candidate_id in attrs.get(candidate_attribute, ()):
                if candidate_id in fasta_ids:
                    matches[(row_feature, candidate_attribute)].add(candidate_id)

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

    selected_feature, selected_attribute = selected_pair
    intervals: dict[str, Gene] = {}
    for seqid, row_feature, start, end, strand, attrs in _iter_gff(gff_path):
        if row_feature != selected_feature:
            continue
        for gene_id in attrs.get(selected_attribute, ()):
            if gene_id not in selected_matches:
                continue
            previous = intervals.get(gene_id)
            if previous is None:
                intervals[gene_id] = Gene(seqid, start, end, gene_id, strand)
            else:
                if previous.seqid != seqid:
                    raise ValueError(f"GFF identifier {gene_id!r} occurs on multiple sequences")
                merged_strand = previous.strand if previous.strand == strand else "."
                intervals[gene_id] = Gene(
                    seqid,
                    min(previous.start, start),
                    max(previous.end, end),
                    gene_id,
                    merged_strand,
                )

    genes = tuple(sorted(intervals.values(), key=lambda gene: (natural_key(gene.seqid), gene.start, gene.end, gene.gene_id)))
    return AnnotationMapping(
        feature=selected_feature,
        attribute=selected_attribute,
        genes=genes,
        fasta_gene_count=len(fasta_ids),
        matched_gene_count=len(genes),
    )


def write_bed(genes: Iterable[Gene], path: str | Path) -> None:
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        handle.writelines(
            f"{gene.seqid}\t{gene.start}\t{gene.end}\t{gene.gene_id}\t0\t{gene.strand}\n" for gene in genes
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
            if not columns[0]:
                raise ValueError(f"Empty BED sequence identifier at {path}:{line_number}")
            if not gene_id:
                raise ValueError(f"Empty BED gene identifier at {path}:{line_number}")
            if gene_id in seen:
                raise ValueError(f"Duplicate BED identifier {gene_id!r} in {path}")
            seen.add(gene_id)
            try:
                start, end = int(columns[1]), int(columns[2])
            except ValueError as exc:
                raise ValueError(f"Invalid BED coordinates at {path}:{line_number}") from exc
            if start < 0 or end <= start:
                raise ValueError(f"Invalid BED interval at {path}:{line_number}")
            strand = columns[5] if len(columns) >= 6 else "."
            genes.append(Gene(columns[0], start, end, gene_id, strand))
    if not genes:
        raise ValueError(f"No BED records found in {path}")
    return tuple(sorted(genes, key=lambda gene: (natural_key(gene.seqid), gene.start, gene.end, gene.gene_id)))


def _first_member(candidates: Iterable[str], identifiers: set[str]) -> str | None:
    return next((candidate for candidate in candidates if candidate in identifiers), None)


def read_jcvi_pairs(path: str | Path, target_ids: set[str], query_ids: set[str]) -> tuple[tuple[str, str], ...]:
    pairs: set[tuple[str, str]] = set()
    unmatched = 0
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 2:
                continue
            left, right = columns[0], columns[1]
            if left in target_ids and right in query_ids:
                pairs.add((left, right))
            elif right in target_ids and left in query_ids:
                pairs.add((right, left))
            else:
                unmatched += 1
    if not pairs:
        raise ValueError(f"No JCVI anchor pairs matched both BED files ({unmatched} unmatched rows)")
    return tuple(sorted(pairs))


def read_synmap_pairs(path: str | Path, target_ids: set[str], query_ids: set[str]) -> tuple[tuple[str, str], ...]:
    pairs: set[tuple[str, str]] = set()
    unmatched = 0
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.replace("||", "\t").rstrip("\n").split("\t")
            if len(columns) < 20:
                unmatched += 1
                continue
            left_candidates = (columns[4], columns[7])
            right_candidates = (columns[16], columns[19])
            left_target = _first_member(left_candidates, target_ids)
            right_query = _first_member(right_candidates, query_ids)
            right_target = _first_member(right_candidates, target_ids)
            left_query = _first_member(left_candidates, query_ids)
            if left_target and right_query:
                pairs.add((left_target, right_query))
            elif right_target and left_query:
                pairs.add((right_target, left_query))
            else:
                unmatched += 1
    if not pairs:
        raise ValueError(f"No SynMap pairs matched both BED files ({unmatched} unmatched rows)")
    return tuple(sorted(pairs))


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
    if format_name == "auto":
        format_name = detect_synteny_format(path)
    if format_name == "jcvi":
        return read_jcvi_pairs(path, target_ids, query_ids)
    if format_name == "synmap":
        return read_synmap_pairs(path, target_ids, query_ids)
    raise ValueError(f"Unsupported synteny format: {format_name}")

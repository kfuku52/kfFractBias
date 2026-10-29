"""Copy-aware, reference-aligned fractionation reports from closed JCVI anchors.

Local copy labels are supplied independently of the measured outcome. Native
segments belonging to one copy are combined by logical OR before windowing.
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import __version__
from .io import Gene, open_text, read_bed, sha256_file
from .provenance import source_metadata
from .run import analysis_run, validate_prefix, validate_separation


def table(path: Path, required: set[str]) -> list[dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"Missing or duplicate table headers: {path}")
        if not required <= set(reader.fieldnames):
            raise ValueError(f"Missing columns {sorted(required - set(reader.fieldnames))}: {path}")
        rows = list(reader)
    if not rows or any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError(f"Empty or ragged table: {path}")
    return rows


@dataclass(frozen=True)
class ComparisonInput:
    analysis_id: str
    species: str
    expected_copies: int
    target_bed: Path
    query_bed: Path
    synteny_file: Path
    mapping_file: Path | None


def read_manifest(path: Path) -> list[ComparisonInput]:
    required = {
        "analysis_id",
        "species",
        "expected_copies",
        "target_bed",
        "query_bed",
        "synteny_file",
    }
    rows = table(path, required)
    if set(rows[0]) - required - {"mapping_file"}:
        raise ValueError("Unsupported comparison manifest columns")
    result = []
    seen = set()
    for row in rows:
        name = row["analysis_id"]
        validate_prefix(name)
        if name in seen or not row["species"].strip():
            raise ValueError("Duplicate analysis ID or empty species")
        seen.add(name)
        copies = int(row["expected_copies"])
        if copies < 1 or str(copies) != row["expected_copies"]:
            raise ValueError("expected_copies must be a positive integer")
        paths = {}
        for key in ["target_bed", "query_bed", "synteny_file", "mapping_file"]:
            value = row.get(key, "")
            paths[key] = (path.parent / value).resolve() if value else None
            if key != "mapping_file" and paths[key] is None:
                raise ValueError(f"Missing {key}")
        target, query, synteny = (paths[k] for k in ["target_bed", "query_bed", "synteny_file"])
        assert target is not None and query is not None and synteny is not None
        result.append(
            ComparisonInput(
                name, row["species"], copies, target, query, synteny, paths["mapping_file"]
            )
        )
    return result


def anchor_blocks(path: Path, target: set[str], query: set[str]) -> Iterator[list[tuple[str, str]]]:
    if target & query:
        raise ValueError("Target and query gene identifiers must be disjoint")
    block: list[tuple[str, str]] = []
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("#"):
                if block:
                    yield block
                block = []
            elif line.strip():
                fields = line.split()
                if len(fields) < 2:
                    raise ValueError(f"Invalid JCVI anchor: {path}:{line_number}")
                a, b = fields[:2]
                if a in query and b in target:
                    a, b = b, a
                if a not in target or b not in query:
                    raise ValueError(f"Unknown anchor identifiers: {path}:{line_number}")
                block.append((a, b))
    if block:
        yield block


def build_loci(
    item: ComparisonInput, *, minimum_anchors: int = 10, trim: int = 5, raw: bool = False
) -> tuple[tuple[Gene, ...], list[dict[str, Any]], dict[str, int]]:
    target = read_bed(item.target_bed)
    query = {g.gene_id: g for g in read_bed(item.query_bed)}
    ranks = {}
    chromosomes: dict[str, list[str]] = defaultdict(list)
    for locus in target:
        ranks[locus.gene_id] = (locus.seqid, len(chromosomes[locus.seqid]))
        chromosomes[locus.seqid].append(locus.gene_id)
    mapping = {}
    labels: dict[str, set[str]] = defaultdict(set)
    if item.mapping_file:
        for row in table(item.mapping_file, {"gene_id", "group_id", "copy", "segment_id"}):
            gene, group, copy, segment = (
                row[k] for k in ["gene_id", "group_id", "copy", "segment_id"]
            )
            if gene not in query or gene in mapping or not segment:
                raise ValueError("Unknown/duplicate mapping gene or empty segment")
            if any(not s or any(c.isspace() or c in ":;|" for c in s) for s in [group, copy]):
                raise ValueError("Group and copy labels must be nonempty, without whitespace/:;|")
            mapping[gene] = (group, copy, segment)
            labels[group].add(copy)
        if any(len(v) != item.expected_copies for v in labels.values()):
            raise ValueError("Every mapped group must contain expected_copies distinct labels")
    elif raw:
        mapping = {g: ("all", v.seqid, v.seqid) for g, v in query.items()}
    else:
        raise ValueError("A copy mapping is required for callable profiles")
    called: dict[str, Any] = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    detected: dict[str, Any] = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    anchors = unassigned = interiors = 0
    for block in anchor_blocks(item.synteny_file, set(ranks), set(query)):
        buckets: dict[tuple[str, str, str, str], set[int]] = defaultdict(set)
        for ref, gene in set(block):
            anchors += 1
            if gene not in mapping:
                unassigned += 1
                continue
            group, copy, segment = mapping[gene]
            # Segment identifiers are scoped to a native query sequence.
            segment = query[gene].seqid + ":" + segment
            chrom, rank = ranks[ref]
            detected[ref][group][copy].add(segment)
            buckets[group, copy, chrom, segment].add(rank)
        for (group, copy, chrom, segment), rr in buckets.items():
            if len(rr) < minimum_anchors:
                continue
            lo, hi = min(rr) + trim, max(rr) - trim
            if lo > hi:
                continue
            interiors += 1
            for ref in chromosomes[chrom][lo : hi + 1]:
                called[ref][group][copy].add(segment)
    if not anchors:
        raise ValueError("No JCVI anchors")
    loci: list[dict[str, Any]] = []
    for locus in target:
        ref = locus.gene_id
        candidates = [g for g, copies in called[ref].items() if set(copies) == labels.get(g)]
        group = candidates[0] if len(candidates) == 1 else ""
        calls = {
            f"{g}:{c}": sorted(s) for g, copies in detected[ref].items() for c, s in copies.items()
        }
        if raw:
            values = dict.fromkeys(calls, 1)
            segments = calls
        elif group:
            values = {
                f"{group}:{c}": int(c in detected[ref].get(group, {}))
                for c in sorted(labels[group])
            }
            segments = {f"{group}:{c}": sorted(s) for c, s in called[ref][group].items()}
        else:
            values, segments = {}, {}
        loci.append(
            {
                "target_seqid": locus.seqid,
                "target_gene": ref,
                "target_rank": ranks[ref][1] + 1,
                "group_id": group,
                "callable": int(raw or bool(group)),
                "retained_copies": sum(values.values()),
                "expected_copies": item.expected_copies,
                "copies": values,
                "segments": segments,
                "ambiguous_groups": len(candidates) if len(candidates) > 1 else 0,
            }
        )
    audit = {
        "anchors": anchors,
        "unassigned_anchors": unassigned,
        "callable_interiors": interiors,
        "mapped_query_genes": len(mapping),
        "query_genes": len(query),
        "ambiguous_loci": sum(bool(r["ambiguous_groups"]) for r in loci),
        "above_expected_copy_loci": sum(r["retained_copies"] > item.expected_copies for r in loci),
    }
    return target, loci, audit


def aligned_windows(
    profiles: list[list[dict[str, Any]]],
    *,
    window_size: int,
    step_size: int,
    denominator: str,
    minimum_callable_fraction: float,
) -> tuple[list[list[dict[str, Any]]], list[bool]]:
    if window_size < 1 or step_size < 1 or denominator not in {"all", "syntenic", "callable"}:
        raise ValueError("Invalid window/step/denominator")
    if not 0 < minimum_callable_fraction <= 1:
        raise ValueError("minimum callable fraction must be in (0,1]")

    def identity(rows: list[dict[str, Any]]) -> list[tuple[Any, Any, Any]]:
        return [(r["target_seqid"], r["target_gene"], r["target_rank"]) for r in rows]

    if not profiles or any(identity(p) != identity(profiles[0]) for p in profiles[1:]):
        raise ValueError("Profiles must have identical reference loci, coordinates and order")
    shared = [all(p[i]["callable"] for p in profiles) for i in range(len(profiles[0]))]
    selected = [
        i
        for i in range(len(shared))
        if denominator == "all"
        or (denominator == "callable" and shared[i])
        or (denominator == "syntenic" and all(p[i]["retained_copies"] > 0 for p in profiles))
    ]
    by_chr: dict[str, list[int]] = defaultdict(list)
    for i in selected:
        by_chr[profiles[0][i]["target_seqid"]].append(i)
    results: list[list[dict[str, Any]]] = [[] for _ in profiles]
    for chrom, indices in by_chr.items():
        common = [0]
        totals = [[0] for _ in profiles]
        for i in indices:
            common.append(common[-1] + int(shared[i]))
            for a, p in enumerate(profiles):
                totals[a].append(totals[a][-1] + (p[i]["retained_copies"] if shared[i] else 0))
        for start in range(0, len(indices) - window_size + 1, step_size):
            end = start + window_size
            members = indices[start:end]
            valid = [i for i in members if shared[i]]
            n = common[end] - common[start]
            estimable = n >= math.ceil(window_size * minimum_callable_fraction)
            for a, p in enumerate(profiles):
                count = totals[a][end] - totals[a][start]
                copies: dict[str, float] = {}
                segments: dict[str, str] = {}
                groups = {p[i]["group_id"] for i in valid}
                if estimable and len(groups) == 1:
                    keys = set().union(*(p[i]["copies"] for i in valid))
                    for key in sorted(keys):
                        native = {s for i in valid for s in p[i]["segments"].get(key, [])}
                        if len(native) == 1:
                            copies[key] = 100 * sum(p[i]["copies"].get(key, 0) for i in valid) / n
                            segments[key] = next(iter(native))
                results[a].append(
                    {
                        "target_seqid": chrom,
                        "start_gene": p[members[0]]["target_gene"],
                        "end_gene": p[members[-1]]["target_gene"],
                        "start_rank": p[members[0]]["target_rank"],
                        "end_rank": p[members[-1]]["target_rank"],
                        "window_size": window_size,
                        "common_callable_genes": n,
                        "retained_copies": count if estimable else None,
                        "total_percent": 100 * count / n if estimable else None,
                        "normalised_percent": 100 * count / (n * p[0]["expected_copies"])
                        if estimable
                        else None,
                        "copy_percent": copies,
                        "segment_ids": segments,
                        "status": "estimated" if estimable else "insufficient_common_callable",
                    }
                )
    return results, shared


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("No rows to write")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    k: json.dumps(v, sort_keys=True) if isinstance(v, (dict, list)) else v
                    for k, v in row.items()
                }
            )


def comparison_report(
    manifest: Path,
    output_dir: Path,
    *,
    prefix: str = "comparison",
    window_size: int = 100,
    step_size: int = 1,
    denominator: str = "all",
    minimum_callable_fraction: float = 1.0,
    minimum_anchors: int = 10,
    trim: int = 5,
    raw: bool = False,
    target_seqids: tuple[str, ...] = (),
    plot_config: dict[str, Any] | None = None,
    plot_config_file: Path | None = None,
) -> dict[str, Any]:
    if minimum_anchors < 1 or trim < 0:
        raise ValueError("Invalid callable-interior criteria")
    manifest_digest = sha256_file(manifest)
    items = read_manifest(manifest)
    sources = {"manifest": manifest.resolve()}
    config_digest = None
    if plot_config_file is not None:
        if plot_config is not None:
            raise ValueError("Specify a plot configuration object or file, not both")
        config_digest = sha256_file(plot_config_file)
        plot_config = json.loads(plot_config_file.read_text())
        if not isinstance(plot_config, dict):
            raise ValueError("Plot configuration must be an object")
        plot_config = {
            **plot_config,
            "font_files": [
                str((plot_config_file.parent / p).resolve())
                for p in plot_config.get("font_files", [])
            ],
        }
        sources["plot_config"] = plot_config_file.resolve()
    for item in items:
        for key in ["target_bed", "query_bed", "synteny_file", "mapping_file"]:
            value = getattr(item, key)
            if value is not None:
                sources[item.analysis_id + ":" + key] = value
    config = dict(plot_config or {})
    for i, name in enumerate(config.get("font_files", [])):
        sources[f"font:{i}"] = Path(name).resolve()
    validate_prefix(prefix)
    outputs = {
        key: output_dir.resolve() / f"{prefix}.{suffix}"
        for key, suffix in [
            ("loci", "loci.tsv"),
            ("windows", "report_windows.tsv"),
            ("summary", "report.json"),
        ]
    }
    for kind in ["normalised", "raw"]:
        for ext in ["pdf", "png", "svg"]:
            outputs[kind + "_" + ext] = output_dir.resolve() / f"{prefix}.{kind}.{ext}"
    outputs["all_pdf"] = output_dir.resolve() / f"{prefix}.profiles_all.pdf"
    validate_separation(sources, outputs)
    with analysis_run(output_dir, prefix, sources) as run:
        frozen = run.capture_inputs(sources)
        if sha256_file(frozen["manifest"]) != manifest_digest:
            raise RuntimeError("Manifest changed during preparation")
        if config_digest is not None and sha256_file(frozen["plot_config"]) != config_digest:
            raise RuntimeError("Plot configuration changed during preparation")
        profiles, audits, reference = [], [], None
        for item in items:
            values = {
                k: frozen[item.analysis_id + ":" + k]
                for k in ["target_bed", "query_bed", "synteny_file", "mapping_file"]
                if getattr(item, k) is not None
            }
            captured = ComparisonInput(
                item.analysis_id,
                item.species,
                item.expected_copies,
                values["target_bed"],
                values["query_bed"],
                values["synteny_file"],
                values.get("mapping_file"),
            )
            genes, loci, audit = build_loci(
                captured, minimum_anchors=minimum_anchors, trim=trim, raw=raw
            )
            if target_seqids:
                if set(target_seqids) - {g.seqid for g in genes}:
                    raise ValueError("Unknown requested reference sequence")
                loci = [r for r in loci if r["target_seqid"] in target_seqids]
                genes = tuple(g for g in genes if g.seqid in target_seqids)
            if reference is not None and genes != reference:
                raise ValueError("Reference BED coordinates differ across analyses")
            reference = genes
            profiles.append(loci)
            audits.append({"analysis_id": item.analysis_id, **audit})
        windows, shared = aligned_windows(
            profiles,
            window_size=window_size,
            step_size=step_size,
            denominator=denominator,
            minimum_callable_fraction=minimum_callable_fraction,
        )
        if not windows[0]:
            raise ValueError("No complete reference window")
        stage = {outputs[k]: run.staging_dir / outputs[k].name for k in outputs}
        commit_paths: dict[Path, Path | None] = dict(stage)
        write_rows(
            stage[outputs["loci"]],
            [
                {"analysis_id": item.analysis_id, **r, "common_callable": int(shared[i])}
                for item, p in zip(items, profiles, strict=True)
                for i, r in enumerate(p)
            ],
        )
        write_rows(
            stage[outputs["windows"]],
            [
                {"analysis_id": item.analysis_id, **r}
                for item, p in zip(items, windows, strict=True)
                for r in p
            ],
        )
        plot = {}
        if plot_config is not None:
            from .comparison_plot import plot_comparison

            plot = plot_comparison(items, windows, commit_paths, outputs, config)
        else:
            for key in list(outputs):
                if key not in {"loci", "windows", "summary"}:
                    commit_paths[outputs[key]] = None
        summary = {
            "schema_version": 1,
            "program": "kfFractBias",
            "program_version": __version__,
            "analysis_mode": "reference_aligned_copy_profiles",
            "source": source_metadata(),
            "inputs": run.inputs,
            "parameters": {
                "window_size": window_size,
                "step_size": step_size,
                "denominator": denominator,
                "minimum_callable_fraction": minimum_callable_fraction,
                "minimum_anchors": minimum_anchors,
                "trim": trim,
                "raw_unmasked": raw,
                "target_seqids": target_seqids,
                "plot_config": config if plot_config is not None else None,
            },
            "analyses": [
                {
                    "analysis_id": i.analysis_id,
                    "species": i.species,
                    "expected_copies": i.expected_copies,
                }
                for i in items
            ],
            "audits": audits,
            "reference_loci": len(shared),
            "common_callable_loci": sum(shared),
            "windows_per_analysis": len(windows[0]),
            "estimable_windows_per_analysis": sum(r["status"] == "estimated" for r in windows[0]),
            "interpretation": "Exploratory syntelog detection conditional on annotation, alignments and local assignments; not validated ancestral gene loss.",
            "plot": plot,
            "output_sha256": {
                k: sha256_file(v)
                for k, p in outputs.items()
                if (v := commit_paths[p]) is not None and k != "summary"
            },
        }
        stage[outputs["summary"]].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        run.commit(commit_paths)
    return summary

from __future__ import annotations

import csv
import json
import platform
import re
import sys
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, TypeVar

from . import __version__
from .io import Gene, detect_synteny_format, natural_key, parse_synteny_pairs, read_bed, sha256_file
from .profiles import GENE_FIELDS, WINDOW_FIELDS, GeneRow, RetentionProfile, WindowRow
from .provenance import source_metadata
from .run import RunContext, analysis_run, output_paths, validate_prefix, validate_separation


@dataclass(frozen=True)
class AnalysisConfig:
    synteny_path: Path
    synteny_format: str
    target_bed: Path
    query_bed: Path
    output_dir: Path
    prefix: str = "kffractbias"
    target_name: str = "target"
    query_name: str = "query"
    window_size: int = 100
    step_size: int = 1
    denominator: str = "all"
    analysis_mode: str = "pairwise_fractionation_bias"
    target_seqids: tuple[str, ...] = ()
    query_seqids: tuple[str, ...] = ()
    exclude_seqid_regex: str = ""
    include_unmatched_query_seqids: bool = False
    make_plot: bool = True
    collect_rows: bool = True
    keep_failed_work: bool = False
    metadata: dict[str, Any] = field(default_factory=dict)
    additional_inputs: dict[str, Path] = field(default_factory=dict)


@dataclass(frozen=True)
class AnalysisResult:
    genes_path: Path
    windows_path: Path
    summary_path: Path
    pdf_path: Path | None
    png_path: Path | None
    gene_rows: tuple[GeneRow, ...]
    window_rows: tuple[WindowRow, ...]


Row = TypeVar("Row", bound=Mapping[str, object])


def _write_rows(
    path: Path, fieldnames: tuple[str, ...], rows: Iterable[Row], *, collect: bool
) -> tuple[Row, ...]:
    collected: list[Row] = []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fieldnames)
        for row in rows:
            writer.writerow([row[field] for field in fieldnames])
            if collect:
                collected.append(row)
    return tuple(collected)


def _write_json(path: Path, value: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _parse_seqids(values: tuple[str, ...]) -> set[str]:
    parsed: set[str] = set()
    for value in values:
        parsed.update(item.strip() for item in value.split(",") if item.strip())
    return parsed


def _filter_genes(
    genes: tuple[Gene, ...],
    requested_seqids: tuple[str, ...],
    exclude_seqid_regex: str,
    label: str,
) -> tuple[Gene, ...]:
    requested = _parse_seqids(requested_seqids)
    available = {gene.seqid for gene in genes}
    missing = requested - available
    if missing:
        raise ValueError(
            f"Unknown {label} sequence identifiers: {', '.join(sorted(missing, key=natural_key))}"
        )
    try:
        pattern = re.compile(exclude_seqid_regex) if exclude_seqid_regex else None
    except re.error as exc:
        raise ValueError(f"Invalid sequence exclusion regex: {exc}") from exc
    selected = tuple(
        gene
        for gene in genes
        if (not requested or gene.seqid in requested)
        and (pattern is None or pattern.search(gene.seqid) is None)
    )
    if not selected:
        raise ValueError(f"No {label} genes remained after sequence filtering")
    return selected


def _validate_config(config: AnalysisConfig) -> None:
    if config.window_size < 1:
        raise ValueError("window size must be at least 1")
    if config.step_size < 1:
        raise ValueError("step size must be at least 1")
    if config.denominator not in {"all", "syntenic"}:
        raise ValueError("denominator must be 'all' or 'syntenic'")
    if config.analysis_mode not in {"pairwise_fractionation_bias", "self_synteny_retention"}:
        raise ValueError(
            "analysis mode must be 'pairwise_fractionation_bias' or 'self_synteny_retention'"
        )
    validate_prefix(config.prefix)
    reserved_inputs = {"synteny", "target_bed", "query_bed"}
    conflicts = reserved_inputs.intersection(config.additional_inputs)
    if conflicts:
        raise ValueError(f"Additional input labels are reserved: {', '.join(sorted(conflicts))}")


def _input_paths(config: AnalysisConfig) -> dict[str, Path]:
    paths = {
        "synteny": Path(config.synteny_path).resolve(),
        "target_bed": Path(config.target_bed).resolve(),
        "query_bed": Path(config.query_bed).resolve(),
    }
    paths.update(
        {label: Path(path).resolve() for label, path in sorted(config.additional_inputs.items())}
    )
    return paths


def _runtime_metadata() -> dict[str, Any]:
    packages: dict[str, str | None] = {}
    for distribution in ("kffractbias", "jcvi", "matplotlib", "ortools"):
        try:
            packages[distribution] = version(distribution)
        except PackageNotFoundError:
            packages[distribution] = None
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": packages,
        "source": source_metadata(),
    }


def _validate_pairwise_identifiers(
    analysis_mode: str,
    target_by_id: dict[str, Gene],
    query_by_id: dict[str, Gene],
) -> None:
    if analysis_mode != "pairwise_fractionation_bias":
        return
    overlapping_ids = set(target_by_id).intersection(query_by_id)
    if overlapping_ids:
        examples = ", ".join(sorted(overlapping_ids, key=natural_key)[:10])
        raise ValueError(
            "Pairwise target and query BED gene identifiers must be disjoint; "
            f"found {len(overlapping_ids)} overlapping identifier(s), including: {examples}"
        )


def _select_filtered_pairs(
    parsed_pairs: tuple[tuple[str, str], ...],
    target_by_id: dict[str, Gene],
    query_by_id: dict[str, Gene],
) -> tuple[tuple[str, str], ...]:
    selected = tuple(
        (target_id, query_id)
        for target_id, query_id in parsed_pairs
        if target_id in target_by_id and query_id in query_by_id
    )
    if not selected:
        raise ValueError("No synteny pairs remained after target/query sequence filtering")
    return selected


def _prepare_directed_pairs(
    analysis_mode: str,
    input_pairs: tuple[tuple[str, str], ...],
    target_by_id: dict[str, Gene],
    query_by_id: dict[str, Gene],
) -> tuple[tuple[tuple[str, str], ...], tuple[tuple[str, str], ...], dict[str, int]]:
    if analysis_mode != "self_synteny_retention":
        return input_pairs, input_pairs, {}

    target_coordinates = {
        gene_id: (gene.seqid, gene.start, gene.end, gene.strand)
        for gene_id, gene in target_by_id.items()
    }
    query_coordinates = {
        gene_id: (gene.seqid, gene.start, gene.end, gene.strand)
        for gene_id, gene in query_by_id.items()
    }
    if target_coordinates != query_coordinates:
        raise ValueError(
            "Self-synteny retention requires identical target and query BED gene sets and coordinates"
        )
    identity_pair_count = sum(target_id == query_id for target_id, query_id in input_pairs)
    canonical_pairs = {
        (target_id, query_id) if target_id < query_id else (query_id, target_id)
        for target_id, query_id in input_pairs
        if target_id != query_id
    }
    if not canonical_pairs:
        raise ValueError("No non-identity self-synteny pairs remained after filtering")
    nonidentity_pair_count = len(input_pairs) - identity_pair_count
    pairs = tuple(sorted(canonical_pairs))
    directed_pairs = tuple(
        sorted(
            (source, destination)
            for left, right in pairs
            for source, destination in ((left, right), (right, left))
        )
    )
    pair_counts = {
        "input_synteny_pair_count": len(input_pairs),
        "removed_identity_pair_count": identity_pair_count,
        "removed_mirrored_pair_count": nonidentity_pair_count - len(canonical_pairs),
        "directed_synteny_pair_count": len(directed_pairs),
        "intrachromosomal_pair_count": sum(
            target_by_id[left].seqid == target_by_id[right].seqid for left, right in pairs
        ),
        "interchromosomal_pair_count": sum(
            target_by_id[left].seqid != target_by_id[right].seqid for left, right in pairs
        ),
    }
    return pairs, directed_pairs, pair_counts


def preflight_analysis(
    config: AnalysisConfig, target_genes: tuple[Gene, ...], query_genes: tuple[Gene, ...]
) -> dict[str, int]:
    """Validate selectors and estimate dense output before expensive alignment."""
    _validate_config(config)
    target = _filter_genes(target_genes, config.target_seqids, config.exclude_seqid_regex, "target")
    query = _filter_genes(query_genes, config.query_seqids, config.exclude_seqid_regex, "query")
    _validate_pairwise_identifiers(
        config.analysis_mode,
        {g.gene_id: g for g in target_genes},
        {g.gene_id: g for g in query_genes},
    )
    grouped: dict[str, list[Gene]] = defaultdict(list)
    for gene in target:
        grouped[gene.seqid].append(gene)
    profile = RetentionProfile(
        grouped,
        sorted({g.seqid for g in query}, key=natural_key),
        {},
        config.window_size,
        config.step_size,
    )
    genes, windows = profile.row_counts()
    return {"gene_table_rows_upper_bound": genes, "window_table_rows_upper_bound": windows}


def calculate_fractionation_bias(
    config: AnalysisConfig, *, run: RunContext | None = None
) -> AnalysisResult:
    _validate_config(config)
    inputs = _input_paths(config)
    outputs = output_paths(config.output_dir, config.prefix)
    validate_separation(inputs, outputs)
    if run is None:
        with analysis_run(
            config.output_dir, config.prefix, inputs, keep_failed_work=config.keep_failed_work
        ) as owned_run:
            return calculate_fractionation_bias(config, run=owned_run)
    run.require_active(config.output_dir, config.prefix)
    inputs = run.capture_inputs(inputs)
    if config.make_plot:
        from .plotting import preflight_plot

        preflight_plot()
    return _calculate(config, run, inputs, outputs)


def _calculate(
    config: AnalysisConfig, run: RunContext, inputs: dict[str, Path], outputs: dict[str, Path]
) -> AnalysisResult:
    metadata = dict(config.metadata)
    if config.analysis_mode == "self_synteny_retention":
        metadata.setdefault(
            "interpretation",
            "Within-genome self-synteny retention conditional on extant annotated genes; not an outgroup-based fractionation-bias estimate.",
        )
    with run.stage("parse_and_prepare_profile"):
        all_target_genes = read_bed(inputs["target_bed"])
        all_query_genes = read_bed(inputs["query_bed"])
        target_genes = _filter_genes(
            all_target_genes,
            config.target_seqids,
            config.exclude_seqid_regex,
            "target",
        )
        query_genes = _filter_genes(
            all_query_genes,
            config.query_seqids,
            config.exclude_seqid_regex,
            "query",
        )
        all_target_by_id = {gene.gene_id: gene for gene in all_target_genes}
        all_query_by_id = {gene.gene_id: gene for gene in all_query_genes}
        target_by_id = {gene.gene_id: gene for gene in target_genes}
        query_by_id = {gene.gene_id: gene for gene in query_genes}
        _validate_pairwise_identifiers(config.analysis_mode, all_target_by_id, all_query_by_id)

        resolved_synteny_format = (
            detect_synteny_format(inputs["synteny"])
            if config.synteny_format == "auto"
            else config.synteny_format
        )
        parsed_synteny = parse_synteny_pairs(
            inputs["synteny"],
            resolved_synteny_format,
            set(all_target_by_id),
            set(all_query_by_id),
            allow_ambiguous_orientation=config.analysis_mode == "self_synteny_retention",
        )
        input_pairs = _select_filtered_pairs(
            parsed_synteny.pairs,
            target_by_id,
            query_by_id,
        )
        pairs, directed_pairs, pair_counts = _prepare_directed_pairs(
            config.analysis_mode,
            input_pairs,
            target_by_id,
            query_by_id,
        )

        mappings: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
        for target_id, query_id in directed_pairs:
            mappings[target_id][query_by_id[query_id].seqid].add(query_id)

        represented_query_seqids = {query_by_id[query_id].seqid for _, query_id in directed_pairs}
        include_all_query_seqids = config.include_unmatched_query_seqids or bool(
            config.query_seqids
        )
        query_seqids = (
            {gene.seqid for gene in query_genes}
            if include_all_query_seqids
            else represented_query_seqids
        )
        ordered_query_seqids = sorted(query_seqids, key=natural_key)
        target_by_seqid: dict[str, list[Gene]] = defaultdict(list)
        for gene in target_genes:
            if config.denominator == "syntenic" and not any(
                mappings[gene.gene_id].get(query_seqid) for query_seqid in ordered_query_seqids
            ):
                continue
            target_by_seqid[gene.seqid].append(gene)

        profile = RetentionProfile(
            target_by_seqid, ordered_query_seqids, mappings, config.window_size, config.step_size
        )
    gene_count, window_count = profile.row_counts()
    staging_dir = run.staging_dir
    staged_genes = staging_dir / outputs["genes"].name
    staged_windows = staging_dir / outputs["windows"].name
    staged_summary = staging_dir / outputs["summary"].name
    staged_pdf = staging_dir / outputs["plot_pdf"].name
    staged_png = staging_dir / outputs["plot_png"].name
    with run.stage("tables"):
        gene_rows = _write_rows(
            staged_genes, GENE_FIELDS, profile.gene_rows(), collect=config.collect_rows
        )
        window_rows = _write_rows(
            staged_windows, WINDOW_FIELDS, profile.window_rows(), collect=config.collect_rows
        )

    pdf_path: Path | None = None
    png_path: Path | None = None
    plot_metadata: dict[str, int] = {}
    if config.make_plot:
        from .plotting import plot_windows

        with run.stage("plot"):
            plot_metadata = plot_windows(
                profile.window_rows(),
                staged_pdf,
                staged_png,
                target_name=config.target_name,
                query_name=config.query_name,
                window_size=config.window_size,
            )
        pdf_path = outputs["plot_pdf"]
        png_path = outputs["plot_png"]
    run.verify_inputs()
    output_hashes = {
        "genes": sha256_file(staged_genes),
        "windows": sha256_file(staged_windows),
        "plot_pdf": sha256_file(staged_pdf) if config.make_plot else None,
        "plot_png": sha256_file(staged_png) if config.make_plot else None,
    }
    summary = {
        "schema_version": 3,
        "program": "kfFractBias",
        "program_version": __version__,
        "analysis_mode": config.analysis_mode,
        "target_name": config.target_name,
        "query_name": config.query_name,
        "runtime": _runtime_metadata(),
        "parameters": {
            "window_size": config.window_size,
            "step_size": config.step_size,
            "denominator": config.denominator,
            "synteny_format": resolved_synteny_format,
            "requested_synteny_format": config.synteny_format,
            "target_seqids": sorted({gene.seqid for gene in target_genes}, key=natural_key),
            "query_seqids": ordered_query_seqids,
            "exclude_seqid_regex": config.exclude_seqid_regex,
            "include_unmatched_query_seqids": config.include_unmatched_query_seqids,
        },
        "counts": {
            "input_target_gene_count": len(all_target_genes),
            "input_query_gene_count": len(all_query_genes),
            "target_gene_count": len(target_genes),
            "analyzed_target_gene_count": sum(map(len, target_by_seqid.values())),
            "query_gene_count": len(query_genes),
            "analyzed_query_sequence_count": len(ordered_query_seqids),
            "input_synteny_record_count": parsed_synteny.record_count,
            "duplicate_synteny_pair_count": parsed_synteny.duplicate_pair_count,
            "sequence_filtered_synteny_pair_count": len(parsed_synteny.pairs) - len(input_pairs),
            "synteny_pair_count": len(pairs),
            "gene_table_row_count": gene_count,
            "window_table_row_count": window_count,
            **pair_counts,
        },
        "inputs": run.inputs,
        "outputs": {
            "genes": str(outputs["genes"]),
            "windows": str(outputs["windows"]),
            "plot_pdf": str(pdf_path) if pdf_path else None,
            "plot_png": str(png_path) if png_path else None,
        },
        "output_sha256": output_hashes,
        "timings_seconds": run.timings,
        "plot": plot_metadata,
        "metadata": run.published_metadata(metadata),
    }
    _write_json(staged_summary, summary)
    run.commit(
        {
            outputs["genes"]: staged_genes,
            outputs["windows"]: staged_windows,
            outputs["plot_pdf"]: staged_pdf if config.make_plot else None,
            outputs["plot_png"]: staged_png if config.make_plot else None,
            outputs["summary"]: staged_summary,
        }
    )
    return AnalysisResult(
        genes_path=outputs["genes"],
        windows_path=outputs["windows"],
        summary_path=outputs["summary"],
        pdf_path=pdf_path,
        png_path=png_path,
        gene_rows=gene_rows,
        window_rows=window_rows,
    )

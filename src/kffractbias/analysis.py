from __future__ import annotations

import csv
import json
import os
import re
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from . import __version__
from .io import Gene, detect_synteny_format, natural_key, read_bed, read_synteny_pairs, sha256_file


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
    target_seqids: tuple[str, ...] = ()
    query_seqids: tuple[str, ...] = ()
    exclude_seqid_regex: str = ""
    make_plot: bool = True
    metadata: dict[str, Any] = field(default_factory=dict)
    additional_inputs: dict[str, Path] = field(default_factory=dict)


@dataclass(frozen=True)
class AnalysisResult:
    genes_path: Path
    windows_path: Path
    summary_path: Path
    pdf_path: Path | None
    png_path: Path | None
    gene_rows: tuple[dict[str, Any], ...]
    window_rows: tuple[dict[str, Any], ...]


def _atomic_write_rows(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def _atomic_write_json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


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
        raise ValueError(f"Unknown {label} sequence identifiers: {', '.join(sorted(missing, key=natural_key))}")
    try:
        pattern = re.compile(exclude_seqid_regex) if exclude_seqid_regex else None
    except re.error as exc:
        raise ValueError(f"Invalid sequence exclusion regex: {exc}") from exc
    selected = tuple(
        gene
        for gene in genes
        if (not requested or gene.seqid in requested) and (pattern is None or pattern.search(gene.seqid) is None)
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
    if (
        not config.prefix
        or Path(config.prefix).name != config.prefix
        or config.prefix in {".", ".."}
        or any(ord(character) < 32 or ord(character) == 127 for character in config.prefix)
    ):
        raise ValueError("prefix must be a non-empty filename component")
    reserved_inputs = {"synteny", "target_bed", "query_bed"}
    conflicts = reserved_inputs.intersection(config.additional_inputs)
    if conflicts:
        raise ValueError(f"Additional input labels are reserved: {', '.join(sorted(conflicts))}")


def calculate_fractionation_bias(config: AnalysisConfig) -> AnalysisResult:
    _validate_config(config)
    target_genes = _filter_genes(
        read_bed(config.target_bed), config.target_seqids, config.exclude_seqid_regex, "target"
    )
    query_genes = _filter_genes(
        read_bed(config.query_bed), config.query_seqids, config.exclude_seqid_regex, "query"
    )
    target_by_id = {gene.gene_id: gene for gene in target_genes}
    query_by_id = {gene.gene_id: gene for gene in query_genes}
    resolved_synteny_format = (
        detect_synteny_format(config.synteny_path) if config.synteny_format == "auto" else config.synteny_format
    )
    pairs = read_synteny_pairs(
        config.synteny_path,
        resolved_synteny_format,
        set(target_by_id),
        set(query_by_id),
    )

    mappings: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for target_id, query_id in pairs:
        mappings[target_id][query_by_id[query_id].seqid].add(query_id)

    ordered_query_seqids = sorted({gene.seqid for gene in query_genes}, key=natural_key)
    target_by_seqid: dict[str, list[Gene]] = defaultdict(list)
    for gene in target_genes:
        if config.denominator == "syntenic" and not any(
            mappings[gene.gene_id].get(query_seqid) for query_seqid in ordered_query_seqids
        ):
            continue
        target_by_seqid[gene.seqid].append(gene)

    gene_rows: list[dict[str, Any]] = []
    window_rows: list[dict[str, Any]] = []
    for target_seqid in sorted(target_by_seqid, key=natural_key):
        genes = target_by_seqid[target_seqid]
        for rank, gene in enumerate(genes, start=1):
            for query_seqid in ordered_query_seqids:
                query_gene_ids = sorted(mappings[gene.gene_id].get(query_seqid, ()), key=natural_key)
                gene_rows.append(
                    {
                        "target_seqid": target_seqid,
                        "target_gene": gene.gene_id,
                        "target_rank": rank,
                        "query_seqid": query_seqid,
                        "retained": int(bool(query_gene_ids)),
                        "query_genes": ";".join(query_gene_ids),
                    }
                )

        if len(genes) < config.window_size:
            continue
        for query_seqid in ordered_query_seqids:
            presence = [int(bool(mappings[gene.gene_id].get(query_seqid))) for gene in genes]
            for start_index in range(0, len(genes) - config.window_size + 1, config.step_size):
                end_index = start_index + config.window_size
                retained_count = sum(presence[start_index:end_index])
                fraction = retained_count / config.window_size
                window_rows.append(
                    {
                        "target_seqid": target_seqid,
                        "query_seqid": query_seqid,
                        "window_index": start_index // config.step_size + 1,
                        "start_rank": start_index + 1,
                        "end_rank": end_index,
                        "start_gene": genes[start_index].gene_id,
                        "end_gene": genes[end_index - 1].gene_id,
                        "retained_count": retained_count,
                        "window_size": config.window_size,
                        "retention_fraction": f"{fraction:.10g}",
                        "retention_percent": f"{fraction * 100:.10g}",
                    }
                )

    config.output_dir.mkdir(parents=True, exist_ok=True)
    genes_path = config.output_dir / f"{config.prefix}.genes.tsv"
    windows_path = config.output_dir / f"{config.prefix}.windows.tsv"
    summary_path = config.output_dir / f"{config.prefix}.summary.json"
    _atomic_write_rows(
        genes_path,
        ["target_seqid", "target_gene", "target_rank", "query_seqid", "retained", "query_genes"],
        gene_rows,
    )
    _atomic_write_rows(
        windows_path,
        [
            "target_seqid",
            "query_seqid",
            "window_index",
            "start_rank",
            "end_rank",
            "start_gene",
            "end_gene",
            "retained_count",
            "window_size",
            "retention_fraction",
            "retention_percent",
        ],
        window_rows,
    )

    pdf_path: Path | None = None
    png_path: Path | None = None
    if config.make_plot:
        from .plotting import plot_windows

        pdf_path = config.output_dir / f"{config.prefix}.plot.pdf"
        png_path = config.output_dir / f"{config.prefix}.plot.png"
        plot_windows(
            window_rows,
            pdf_path,
            png_path,
            target_name=config.target_name,
            query_name=config.query_name,
            window_size=config.window_size,
        )

    inputs = {
        "synteny": {"path": str(config.synteny_path), "sha256": sha256_file(config.synteny_path)},
        "target_bed": {"path": str(config.target_bed), "sha256": sha256_file(config.target_bed)},
        "query_bed": {"path": str(config.query_bed), "sha256": sha256_file(config.query_bed)},
    }
    inputs.update(
        {
            label: {"path": str(path), "sha256": sha256_file(path)}
            for label, path in sorted(config.additional_inputs.items())
        }
    )
    summary = {
        "schema_version": 1,
        "program": "kfFractBias",
        "program_version": __version__,
        "target_name": config.target_name,
        "query_name": config.query_name,
        "parameters": {
            "window_size": config.window_size,
            "step_size": config.step_size,
            "denominator": config.denominator,
            "synteny_format": resolved_synteny_format,
            "requested_synteny_format": config.synteny_format,
            "target_seqids": sorted({gene.seqid for gene in target_genes}, key=natural_key),
            "query_seqids": ordered_query_seqids,
            "exclude_seqid_regex": config.exclude_seqid_regex,
        },
        "counts": {
            "target_gene_count": len(target_genes),
            "analyzed_target_gene_count": sum(len(genes) for genes in target_by_seqid.values()),
            "query_gene_count": len(query_genes),
            "synteny_pair_count": len(pairs),
            "gene_table_row_count": len(gene_rows),
            "window_table_row_count": len(window_rows),
        },
        "inputs": inputs,
        "outputs": {
            "genes": str(genes_path),
            "windows": str(windows_path),
            "plot_pdf": str(pdf_path) if pdf_path else None,
            "plot_png": str(png_path) if png_path else None,
        },
        "metadata": config.metadata,
    }
    _atomic_write_json(summary_path, summary)
    return AnalysisResult(
        genes_path=genes_path,
        windows_path=windows_path,
        summary_path=summary_path,
        pdf_path=pdf_path,
        png_path=png_path,
        gene_rows=tuple(gene_rows),
        window_rows=tuple(window_rows),
    )

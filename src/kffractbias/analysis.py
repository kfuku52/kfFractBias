from __future__ import annotations

import csv
import fcntl
import json
import os
import platform
import re
import sys
import tempfile
from collections import defaultdict
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from . import __version__
from .io import Gene, detect_synteny_format, natural_key, parse_synteny_pairs, read_bed, sha256_file


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


def _write_rows(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


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


def _output_paths(config: AnalysisConfig) -> dict[str, Path]:
    output_dir = Path(config.output_dir).resolve()
    return {
        "genes": output_dir / f"{config.prefix}.genes.tsv",
        "windows": output_dir / f"{config.prefix}.windows.tsv",
        "summary": output_dir / f"{config.prefix}.summary.json",
        "plot_pdf": output_dir / f"{config.prefix}.plot.pdf",
        "plot_png": output_dir / f"{config.prefix}.plot.png",
        "lock": output_dir / f".{config.prefix}.lock",
    }


def _paths_refer_to_same_file(left: Path, right: Path) -> bool:
    if left.resolve() == right.resolve():
        return True
    try:
        return left.samefile(right)
    except (FileNotFoundError, OSError):
        return False


def _validate_input_output_separation(inputs: dict[str, Path], outputs: dict[str, Path]) -> None:
    for input_label, input_path in inputs.items():
        for output_label, output_path in outputs.items():
            if _paths_refer_to_same_file(input_path, output_path):
                raise ValueError(
                    f"Input {input_label!r} and output {output_label!r} refer to the same path: {input_path}"
                )


def _snapshot_inputs(inputs: dict[str, Path]) -> dict[str, dict[str, str]]:
    return {
        label: {"path": str(path), "sha256": sha256_file(path)} for label, path in inputs.items()
    }


def _verify_input_snapshots(inputs: dict[str, Path], snapshots: dict[str, dict[str, str]]) -> None:
    changed: list[str] = []
    for label, path in inputs.items():
        try:
            current_hash = sha256_file(path)
        except OSError:
            changed.append(label)
            continue
        if current_hash != snapshots[label]["sha256"]:
            changed.append(label)
    if changed:
        raise RuntimeError(f"Input files changed during analysis: {', '.join(changed)}")


@contextmanager
def _exclusive_output_lock(lock_path: Path) -> Iterator[None]:
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+", encoding="utf-8") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(
                f"Another analysis is writing the same output prefix: {lock_path}"
            ) from exc
        handle.seek(0)
        handle.truncate()
        handle.write(f"pid={os.getpid()}\n")
        handle.flush()
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _remove_file(path: Path) -> None:
    if path.exists() or path.is_symlink():
        path.unlink()


def _commit_outputs(staged: dict[Path, Path | None], staging_dir: Path) -> None:
    backups: dict[Path, Path] = {}
    installed: list[Path] = []
    for final_path in staged:
        if final_path.exists() and final_path.is_dir():
            raise ValueError(f"Output path is an existing directory: {final_path}")
    try:
        for index, (final_path, staged_path) in enumerate(staged.items()):
            if final_path.exists() or final_path.is_symlink():
                backup_path = staging_dir / f"backup-{index}-{final_path.name}"
                os.replace(final_path, backup_path)
                backups[final_path] = backup_path
            if staged_path is not None:
                os.replace(staged_path, final_path)
                installed.append(final_path)
    except BaseException:
        for final_path in installed:
            _remove_file(final_path)
        for final_path, backup_path in backups.items():
            os.replace(backup_path, final_path)
        raise


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
    }


def _build_rows(
    target_by_seqid: dict[str, list[Gene]],
    ordered_query_seqids: list[str],
    mappings: dict[str, dict[str, set[str]]],
    window_size: int,
    step_size: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    gene_rows: list[dict[str, Any]] = []
    window_rows: list[dict[str, Any]] = []
    for target_seqid in sorted(target_by_seqid, key=natural_key):
        genes = target_by_seqid[target_seqid]
        for rank, gene in enumerate(genes, start=1):
            for query_seqid in ordered_query_seqids:
                query_gene_ids = sorted(
                    mappings[gene.gene_id].get(query_seqid, ()), key=natural_key
                )
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

        if len(genes) < window_size:
            continue
        for query_seqid in ordered_query_seqids:
            cumulative = [0]
            for gene in genes:
                cumulative.append(
                    cumulative[-1] + int(bool(mappings[gene.gene_id].get(query_seqid)))
                )
            for start_index in range(0, len(genes) - window_size + 1, step_size):
                end_index = start_index + window_size
                retained_count = cumulative[end_index] - cumulative[start_index]
                fraction = retained_count / window_size
                window_rows.append(
                    {
                        "target_seqid": target_seqid,
                        "query_seqid": query_seqid,
                        "window_index": start_index // step_size + 1,
                        "start_rank": start_index + 1,
                        "end_rank": end_index,
                        "start_gene": genes[start_index].gene_id,
                        "end_gene": genes[end_index - 1].gene_id,
                        "retained_count": retained_count,
                        "window_size": window_size,
                        "retention_fraction": f"{fraction:.10g}",
                        "retention_percent": f"{fraction * 100:.10g}",
                    }
                )
    return gene_rows, window_rows


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


def calculate_fractionation_bias(config: AnalysisConfig) -> AnalysisResult:
    _validate_config(config)
    inputs = _input_paths(config)
    outputs = _output_paths(config)
    _validate_input_output_separation(inputs, outputs)

    with _exclusive_output_lock(outputs["lock"]):
        input_snapshots = _snapshot_inputs(inputs)
        all_target_genes = read_bed(config.target_bed)
        all_query_genes = read_bed(config.query_bed)
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
            detect_synteny_format(config.synteny_path)
            if config.synteny_format == "auto"
            else config.synteny_format
        )
        parsed_synteny = parse_synteny_pairs(
            config.synteny_path,
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

        gene_rows, window_rows = _build_rows(
            target_by_seqid,
            ordered_query_seqids,
            mappings,
            config.window_size,
            config.step_size,
        )

        with tempfile.TemporaryDirectory(
            prefix=f".{config.prefix}.staging-",
            dir=outputs["genes"].parent,
        ) as temporary:
            staging_dir = Path(temporary)
            staged_genes = staging_dir / outputs["genes"].name
            staged_windows = staging_dir / outputs["windows"].name
            staged_summary = staging_dir / outputs["summary"].name
            staged_pdf = staging_dir / outputs["plot_pdf"].name
            staged_png = staging_dir / outputs["plot_png"].name

            _write_rows(
                staged_genes,
                [
                    "target_seqid",
                    "target_gene",
                    "target_rank",
                    "query_seqid",
                    "retained",
                    "query_genes",
                ],
                gene_rows,
            )
            _write_rows(
                staged_windows,
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

                plot_windows(
                    window_rows,
                    staged_pdf,
                    staged_png,
                    target_name=config.target_name,
                    query_name=config.query_name,
                    window_size=config.window_size,
                )
                pdf_path = outputs["plot_pdf"]
                png_path = outputs["plot_png"]

            _verify_input_snapshots(inputs, input_snapshots)
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
                    "analyzed_target_gene_count": sum(
                        len(genes) for genes in target_by_seqid.values()
                    ),
                    "query_gene_count": len(query_genes),
                    "analyzed_query_sequence_count": len(ordered_query_seqids),
                    "input_synteny_record_count": parsed_synteny.record_count,
                    "duplicate_synteny_pair_count": parsed_synteny.duplicate_pair_count,
                    "sequence_filtered_synteny_pair_count": len(parsed_synteny.pairs)
                    - len(input_pairs),
                    "synteny_pair_count": len(pairs),
                    "gene_table_row_count": len(gene_rows),
                    "window_table_row_count": len(window_rows),
                    **pair_counts,
                },
                "inputs": input_snapshots,
                "outputs": {
                    "genes": str(outputs["genes"]),
                    "windows": str(outputs["windows"]),
                    "plot_pdf": str(pdf_path) if pdf_path else None,
                    "plot_png": str(png_path) if png_path else None,
                },
                "output_sha256": output_hashes,
                "metadata": config.metadata,
            }
            _write_json(staged_summary, summary)
            _commit_outputs(
                {
                    outputs["genes"]: staged_genes,
                    outputs["windows"]: staged_windows,
                    outputs["plot_pdf"]: staged_pdf if config.make_plot else None,
                    outputs["plot_png"]: staged_png if config.make_plot else None,
                    outputs["summary"]: staged_summary,
                },
                staging_dir,
            )

    return AnalysisResult(
        genes_path=outputs["genes"],
        windows_path=outputs["windows"],
        summary_path=outputs["summary"],
        pdf_path=outputs["plot_pdf"] if config.make_plot else None,
        png_path=outputs["plot_png"] if config.make_plot else None,
        gene_rows=tuple(gene_rows),
        window_rows=tuple(window_rows),
    )

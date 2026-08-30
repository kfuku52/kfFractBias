from __future__ import annotations

import json
import re
import shutil
import subprocess  # nosec B404
import sys
import time
from collections.abc import Callable
from dataclasses import dataclass, field
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from .io import (
    AnnotationMapping,
    annotation_to_genes,
    iter_fasta,
    select_isoforms,
    write_bed,
)


@dataclass(frozen=True)
class PreparedGenome:
    label: str
    bed_path: Path
    cds_path: Path
    mapping: AnnotationMapping


@dataclass(frozen=True)
class SyntenyRun:
    target: PreparedGenome
    query: PreparedGenome
    anchors_path: Path
    command: tuple[str, ...]
    quota: str
    tool_versions: dict[str, str | None] = field(default_factory=dict)
    commands: tuple[tuple[str, ...], ...] = ()
    blast_task: str | None = None


@dataclass(frozen=True)
class SelfSyntenyRun:
    genome: PreparedGenome
    anchors_path: Path
    commands: tuple[tuple[str, ...], ...]
    depth: int
    quota: str
    self_hit_percent: float
    intrachromosomal_diagonal_bound: int
    tool_versions: dict[str, str | None] = field(default_factory=dict)
    blast_task: str | None = None


def validate_quota(value: str) -> str:
    match = re.fullmatch(r"([1-9]\d*):([1-9]\d*)", value)
    if not match:
        raise ValueError(
            "quota must have the form positive_integer:positive_integer, for example 1:2"
        )
    return value


def _distribution_version(distribution: str) -> str | None:
    try:
        return version(distribution)
    except PackageNotFoundError:
        return None


def _aligner_version(aligner: str) -> str | None:
    executable = "lastal" if aligner == "last" else "blastn"
    resolved = shutil.which(executable)
    if resolved is None:
        return None
    argument = "--version" if aligner == "last" else "-version"
    try:
        output = subprocess.check_output(  # nosec B603
            (resolved, argument),
            stderr=subprocess.STDOUT,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return output.strip().splitlines()[0] if output.strip() else None


def _tool_versions(aligner: str) -> dict[str, str | None]:
    return {
        "jcvi": _distribution_version("jcvi"),
        "ortools": _distribution_version("ortools"),
        aligner: _aligner_version(aligner),
    }


def _run_checked(command: tuple[str, ...], *, cwd: Path, description: str) -> None:
    log_dir = cwd / "logs"
    log_dir.mkdir(exist_ok=True)
    name = f"{len(list(log_dir.glob('*.json'))) + 1:02d}-{re.sub('[^a-z0-9]+', '-', description.lower()).strip('-')}"
    log_path = log_dir / f"{name}.log"
    started = time.perf_counter()
    returncode: int | None = None
    print(f"kffractbias: {description}; log: {log_path}", file=sys.stderr)
    try:
        with log_path.open("wb") as log:
            subprocess.run(command, cwd=cwd, check=True, stdout=log, stderr=subprocess.STDOUT)  # nosec B603
        returncode = 0
    except subprocess.CalledProcessError as exc:
        returncode = exc.returncode
        with log_path.open("rb") as log:
            log.seek(max(0, log_path.stat().st_size - 8192))
            tail = log.read().decode("utf-8", errors="replace").strip()
        raise RuntimeError(
            f"{description} failed with exit code {exc.returncode}\n{tail}\nLog: {log_path}"
        ) from exc
    finally:
        (log_dir / f"{name}.json").write_text(
            json.dumps(
                {
                    "description": description,
                    "command": [part.replace(str(cwd) + "/", "") for part in command],
                    "cwd": ".",
                    "returncode": returncode,
                    "elapsed_seconds": time.perf_counter() - started,
                    "log": str(log_path.relative_to(cwd)),
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )


def preflight_tools(aligner: str) -> None:
    if aligner not in {"last", "blast"}:
        raise ValueError("aligner must be 'last' or 'blast'")
    executables = ("lastal", "lastdb") if aligner == "last" else ("blastn", "makeblastdb")
    missing = [name for name in executables if shutil.which(name) is None]
    if missing:
        raise RuntimeError(f"Missing alignment executables on PATH: {', '.join(missing)}")
    try:
        from jcvi.compara import catalog  # noqa: F401
        from ortools.linear_solver import pywraplp
    except ImportError as exc:
        raise RuntimeError(
            "Synteny generation requires the compare extra: install 'kffractbias[compare]'"
        ) from exc
    if pywraplp.Solver.CreateSolver("SCIP") is None:
        raise RuntimeError("JCVI QUOTA-ALIGN requires the OR-Tools SCIP solver")


def _run_blast(
    target: PreparedGenome, query: PreparedGenome, work_dir: Path, cpus: int, task: str
) -> tuple[tuple[str, ...], ...]:
    if task not in {"blastn", "dc-megablast", "megablast"}:
        raise ValueError("BLAST task must be blastn, dc-megablast, or megablast")
    database = query.cds_path.name
    commands = (
        ("makeblastdb", "-in", database, "-dbtype", "nucl", "-out", database),
        (
            "blastn",
            "-task",
            task,
            "-num_threads",
            str(cpus),
            "-query",
            target.cds_path.name,
            "-db",
            database,
            "-out",
            f"{target.label}.{query.label}.last",
            "-outfmt",
            "6",
            "-max_target_seqs",
            "1000",
            "-evalue",
            "1e-5",
        ),
    )
    for command, description in zip(
        commands, ("BLAST database preparation", "BLAST nucleotide alignment"), strict=True
    ):
        _run_checked(command, cwd=work_dir, description=description)
    return commands


def _write_selected_fasta(source: Path, destination: Path, selected_ids: set[str]) -> None:
    written: set[str] = set()
    with destination.open("w", encoding="utf-8") as output_handle:
        for identifier, header, sequence in iter_fasta(source):
            if identifier in selected_ids:
                written.add(identifier)
                output_handle.write(f">{header}\n")
                for start in range(0, len(sequence), 80):
                    output_handle.write(sequence[start : start + 80] + "\n")
    missing = selected_ids - written
    if missing:
        raise ValueError(f"Failed to copy {len(missing)} selected CDS FASTA records")


def prepare_genome(
    label: str,
    cds_path: str | Path,
    gff_path: str | Path,
    work_dir: str | Path,
    *,
    feature: str | None = None,
    attribute: str | None = None,
    minimum_mapping_fraction: float = 1.0,
    isoform_policy: str = "error",
) -> PreparedGenome:
    work_dir = Path(work_dir).resolve()
    if not 0 < minimum_mapping_fraction <= 1:
        raise ValueError("minimum mapping fraction must be greater than 0 and at most 1")
    lengths = {identifier: len(sequence) for identifier, _header, sequence in iter_fasta(cds_path)}
    fasta_ids = set(lengths)
    mapping = annotation_to_genes(gff_path, fasta_ids, feature=feature, attribute=attribute)
    mapping_fraction = mapping.matched_gene_count / mapping.fasta_gene_count
    if mapping_fraction < minimum_mapping_fraction:
        raise ValueError(
            f"Only {mapping.matched_gene_count}/{mapping.fasta_gene_count} CDS identifiers mapped to {gff_path} "
            f"({mapping_fraction:.1%}); required at least {minimum_mapping_fraction:.1%}"
        )
    mapping = select_isoforms(mapping, lengths, isoform_policy)
    bed_path = work_dir / f"{label}.bed"
    cds_output_path = work_dir / f"{label}.cds"
    write_bed(mapping.genes, bed_path)
    _write_selected_fasta(Path(cds_path), cds_output_path, {gene.gene_id for gene in mapping.genes})
    return PreparedGenome(label, bed_path, cds_output_path, mapping)


def run_pairwise_synteny(
    *,
    target_cds: str | Path,
    target_gff: str | Path,
    query_cds: str | Path,
    query_gff: str | Path,
    work_dir: str | Path,
    quota: str,
    cpus: int,
    cscore: float,
    aligner: str,
    target_feature: str | None = None,
    target_attribute: str | None = None,
    query_feature: str | None = None,
    query_attribute: str | None = None,
    minimum_mapping_fraction: float = 1.0,
    isoform_policy: str = "error",
    before_alignment: Callable[[PreparedGenome, PreparedGenome], None] | None = None,
    blast_task: str = "blastn",
) -> SyntenyRun:
    quota = validate_quota(quota)
    if cpus < 1:
        raise ValueError("cpus must be at least 1")
    if not 0 < cscore <= 1:
        raise ValueError("cscore must be greater than 0 and at most 1")

    work_dir = Path(work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=False)
    target = prepare_genome(
        "target",
        target_cds,
        target_gff,
        work_dir,
        feature=target_feature,
        attribute=target_attribute,
        minimum_mapping_fraction=minimum_mapping_fraction,
        isoform_policy=isoform_policy,
    )
    query = prepare_genome(
        "query",
        query_cds,
        query_gff,
        work_dir,
        feature=query_feature,
        attribute=query_attribute,
        minimum_mapping_fraction=minimum_mapping_fraction,
        isoform_policy=isoform_policy,
    )
    overlapping_ids = {gene.gene_id for gene in target.mapping.genes}.intersection(
        gene.gene_id for gene in query.mapping.genes
    )
    if overlapping_ids:
        examples = ", ".join(sorted(overlapping_ids)[:10])
        raise ValueError(
            "Pairwise target and query CDS/GFF identifiers must be disjoint; "
            f"found {len(overlapping_ids)} overlapping identifier(s), including: {examples}"
        )
    if before_alignment is not None:
        before_alignment(target, query)
    preflight_tools(aligner)
    commands = _run_blast(target, query, work_dir, cpus, blast_task) if aligner == "blast" else ()
    command = (
        sys.executable,
        "-m",
        "jcvi.compara.catalog",
        "ortholog",
        "target",
        "query",
        "--dbtype=nucl",
        "--no_strip_names",
        "--no_dotplot",
        f"--quota={quota}",
        f"--cscore={cscore}",
        f"--cpus={cpus}",
        f"--align_soft={aligner}",
    )
    _run_checked(command, cwd=work_dir, description="JCVI synteny command")
    quota_suffix = quota.replace(":", "x")
    anchors_path = work_dir / f"target.query.lifted.{quota_suffix}.anchors"
    if not anchors_path.is_file() or anchors_path.stat().st_size == 0:
        raise RuntimeError(
            f"JCVI did not create the expected quota-filtered anchors file: {anchors_path}"
        )
    return SyntenyRun(
        target,
        query,
        anchors_path,
        command,
        quota,
        _tool_versions(aligner),
        (*commands, command),
        blast_task if aligner == "blast" else None,
    )


def _validate_self_parameters(
    depth: int,
    cpus: int,
    cscore: float,
    self_hit_percent: float,
    intrachromosomal_diagonal_bound: int,
) -> None:
    if depth < 1:
        raise ValueError("depth must be at least 1")
    if cpus < 1:
        raise ValueError("cpus must be at least 1")
    if not 0 < cscore <= 1:
        raise ValueError("cscore must be greater than 0 and at most 1")
    if not 0 < self_hit_percent <= 100:
        raise ValueError("self hit percent must be greater than 0 and at most 100")
    if intrachromosomal_diagonal_bound < 1:
        raise ValueError("intrachromosomal diagonal bound must be at least 1")


def run_self_synteny(
    *,
    cds: str | Path,
    gff: str | Path,
    work_dir: str | Path,
    depth: int,
    cpus: int,
    cscore: float,
    aligner: str,
    feature: str | None = None,
    attribute: str | None = None,
    minimum_mapping_fraction: float = 1.0,
    self_hit_percent: float = 98.0,
    intrachromosomal_diagonal_bound: int = 300,
    isoform_policy: str = "error",
    before_alignment: Callable[[PreparedGenome, PreparedGenome], None] | None = None,
    blast_task: str = "blastn",
) -> SelfSyntenyRun:
    """Run JCVI alignment with chromosome-aware self chaining and shared-axis quota."""
    _validate_self_parameters(
        depth,
        cpus,
        cscore,
        self_hit_percent,
        intrachromosomal_diagonal_bound,
    )

    work_dir = Path(work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=False)
    genome = prepare_genome(
        "self",
        cds,
        gff,
        work_dir,
        feature=feature,
        attribute=attribute,
        minimum_mapping_fraction=minimum_mapping_fraction,
        isoform_policy=isoform_policy,
    )
    if before_alignment is not None:
        before_alignment(genome, genome)
    preflight_tools(aligner)
    blast_commands = (
        _run_blast(genome, genome, work_dir, cpus, blast_task) if aligner == "blast" else ()
    )
    alignment_command = (
        sys.executable,
        "-m",
        "kffractbias.selfscan",
        "align",
        f"--self-hit-percent={self_hit_percent}",
        f"--cscore={cscore}",
        f"--cpus={cpus}",
        f"--aligner={aligner}",
    )
    _run_checked(alignment_command, cwd=work_dir, description="JCVI self-alignment and filtering")

    filtered_candidates = tuple(work_dir.glob("self.self.last*.inverse.filtered"))
    if len(filtered_candidates) != 1:
        raise RuntimeError("Could not identify the unique JCVI filtered self-alignment")
    filtered_alignment = filtered_candidates[0]
    quota = f"{depth}:{depth}"
    anchors_path = work_dir / f"self.self.lifted.{depth}x{depth}.anchors"
    scan_command = (
        sys.executable,
        "-m",
        "kffractbias.selfscan",
        "scan",
        str(filtered_alignment),
        str(filtered_alignment.with_suffix("")),
        str(genome.bed_path),
        str(anchors_path),
        f"--diagonal-bound={intrachromosomal_diagonal_bound}",
        f"--depth={depth}",
    )
    _run_checked(
        scan_command, cwd=work_dir, description="Chromosome-aware JCVI self scan and QUOTA-ALIGN"
    )
    if not anchors_path.is_file() or anchors_path.stat().st_size == 0:
        raise RuntimeError(
            f"JCVI did not create the expected self quota-filtered anchors file: {anchors_path}"
        )
    commands = [*blast_commands, alignment_command, scan_command]
    return SelfSyntenyRun(
        genome=genome,
        anchors_path=anchors_path,
        commands=tuple(commands),
        depth=depth,
        quota=quota,
        self_hit_percent=self_hit_percent,
        intrachromosomal_diagonal_bound=intrachromosomal_diagonal_bound,
        tool_versions=_tool_versions(aligner),
        blast_task=blast_task if aligner == "blast" else None,
    )

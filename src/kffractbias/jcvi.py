from __future__ import annotations

import re
import shutil
import subprocess  # nosec B404
import sys
from dataclasses import dataclass, field
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from .io import (
    AnnotationMapping,
    annotation_to_genes,
    open_text,
    read_fasta_ids,
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
    try:
        subprocess.run(command, cwd=cwd, check=True)  # nosec B603
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"{description} failed with exit code {exc.returncode}") from exc


def _write_selected_fasta(source: Path, destination: Path, selected_ids: set[str]) -> None:
    written: set[str] = set()
    keep = False
    with (
        open_text(source) as input_handle,
        destination.open("w", encoding="utf-8") as output_handle,
    ):
        for line in input_handle:
            if line.startswith(">"):
                identifier = line[1:].strip().split(None, 1)[0]
                keep = identifier in selected_ids
                if keep:
                    written.add(identifier)
            if keep:
                output_handle.write(line)
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
) -> PreparedGenome:
    work_dir = Path(work_dir)
    if not 0 < minimum_mapping_fraction <= 1:
        raise ValueError("minimum mapping fraction must be greater than 0 and at most 1")
    fasta_ids = read_fasta_ids(cds_path)
    mapping = annotation_to_genes(gff_path, fasta_ids, feature=feature, attribute=attribute)
    mapping_fraction = mapping.matched_gene_count / mapping.fasta_gene_count
    if mapping_fraction < minimum_mapping_fraction:
        raise ValueError(
            f"Only {mapping.matched_gene_count}/{mapping.fasta_gene_count} CDS identifiers mapped to {gff_path} "
            f"({mapping_fraction:.1%}); required at least {minimum_mapping_fraction:.1%}"
        )
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
) -> SyntenyRun:
    quota = validate_quota(quota)
    if cpus < 1:
        raise ValueError("cpus must be at least 1")
    if not 0 < cscore <= 1:
        raise ValueError("cscore must be greater than 0 and at most 1")

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=False)
    target = prepare_genome(
        "target",
        target_cds,
        target_gff,
        work_dir,
        feature=target_feature,
        attribute=target_attribute,
        minimum_mapping_fraction=minimum_mapping_fraction,
    )
    query = prepare_genome(
        "query",
        query_cds,
        query_gff,
        work_dir,
        feature=query_feature,
        attribute=query_attribute,
        minimum_mapping_fraction=minimum_mapping_fraction,
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
    return SyntenyRun(target, query, anchors_path, command, quota, _tool_versions(aligner))


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


def _rerun_self_scan(
    work_dir: Path,
    genome: PreparedGenome,
    intrachromosomal_diagonal_bound: int,
) -> tuple[str, ...]:
    filtered_candidates = tuple(work_dir.glob("self.self.last*.inverse.filtered"))
    if len(filtered_candidates) != 1:
        raise RuntimeError(
            "Could not identify the unique JCVI filtered self-alignment needed to apply "
            f"--diagonal-bound={intrachromosomal_diagonal_bound}"
        )
    filtered_alignment = filtered_candidates[0]
    liftover_alignment = filtered_alignment.with_suffix("")
    anchors = work_dir / "self.self.anchors"
    anchors.unlink(missing_ok=True)
    (work_dir / "self.self.lifted.anchors").unlink(missing_ok=True)
    scan_command = (
        sys.executable,
        "-m",
        "jcvi.compara.synteny",
        "scan",
        str(filtered_alignment),
        str(anchors),
        "--min_size=4",
        "--dist=20",
        f"--liftover={liftover_alignment}",
        f"--intrabound={intrachromosomal_diagonal_bound}",
        "--no_strip_names",
        f"--qbed={genome.bed_path}",
        f"--sbed={genome.bed_path}",
    )
    _run_checked(scan_command, cwd=work_dir, description="JCVI self-synteny scan")
    return scan_command


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
) -> SelfSyntenyRun:
    """Run JCVI in its native self-comparison mode and screen mirrored blocks."""
    _validate_self_parameters(
        depth,
        cpus,
        cscore,
        self_hit_percent,
        intrachromosomal_diagonal_bound,
    )

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=False)
    genome = prepare_genome(
        "self",
        cds,
        gff,
        work_dir,
        feature=feature,
        attribute=attribute,
        minimum_mapping_fraction=minimum_mapping_fraction,
    )
    ortholog_command = (
        sys.executable,
        "-m",
        "jcvi.compara.catalog",
        "ortholog",
        "self",
        "self",
        "--dbtype=nucl",
        "--no_strip_names",
        "--no_dotplot",
        "--ignore_zero_anchor",
        f"--self_remove={self_hit_percent}",
        f"--cscore={cscore}",
        f"--cpus={cpus}",
        f"--align_soft={aligner}",
    )
    _run_checked(ortholog_command, cwd=work_dir, description="JCVI self-synteny command")

    commands: list[tuple[str, ...]] = [ortholog_command]
    lifted_anchors = work_dir / "self.self.lifted.anchors"
    if intrachromosomal_diagonal_bound != 300:
        commands.append(
            _rerun_self_scan(
                work_dir,
                genome,
                intrachromosomal_diagonal_bound,
            )
        )
    if not lifted_anchors.is_file() or lifted_anchors.stat().st_size == 0:
        raise RuntimeError(
            f"JCVI did not create the expected self-synteny anchors file: {lifted_anchors}"
        )

    quota = f"{depth}:{depth}"
    quota_command = (
        sys.executable,
        "-m",
        "jcvi.compara.quota",
        str(lifted_anchors),
        f"--quota={quota}",
        "--self",
        "--screen",
        f"--qbed={genome.bed_path}",
        f"--sbed={genome.bed_path}",
    )
    _run_checked(quota_command, cwd=work_dir, description="JCVI self QUOTA-ALIGN command")

    quota_suffix = quota.replace(":", "x")
    anchors_path = work_dir / f"self.self.lifted.{quota_suffix}.anchors"
    if not anchors_path.is_file() or anchors_path.stat().st_size == 0:
        raise RuntimeError(
            f"JCVI did not create the expected self quota-filtered anchors file: {anchors_path}"
        )
    commands.append(quota_command)
    return SelfSyntenyRun(
        genome=genome,
        anchors_path=anchors_path,
        commands=tuple(commands),
        depth=depth,
        quota=quota,
        self_hit_percent=self_hit_percent,
        intrachromosomal_diagonal_bound=intrachromosomal_diagonal_bound,
        tool_versions=_tool_versions(aligner),
    )

from __future__ import annotations

import argparse
import json
import re
import sys
import tempfile
import time
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

from . import __version__
from .analysis import (
    AnalysisConfig,
    AnalysisResult,
    calculate_fractionation_bias,
    preflight_analysis,
)
from .io import validate_disjoint_identifiers
from .jcvi import (
    PreparedGenome,
    prepare_genome,
    run_pairwise_synteny,
    run_self_synteny,
    validate_quota,
)
from .run import RunContext, analysis_run


def _path(value: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"file not found: {value}")
    return path


def _output_path(value: str) -> Path:
    return Path(value).expanduser().resolve()


def _prefix(value: str) -> str:
    if (
        not value
        or Path(value).name != value
        or value in {".", ".."}
        or any(ord(character) < 32 or ord(character) == 127 for character in value)
    ):
        raise argparse.ArgumentTypeError("prefix must be a non-empty filename component")
    return value


def _positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def _unit_interval(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be greater than 0 and at most 1") from exc
    if not 0 < parsed <= 1:
        raise argparse.ArgumentTypeError("must be greater than 0 and at most 1")
    return parsed


def _percent(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be greater than 0 and at most 100") from exc
    if not 0 < parsed <= 100:
        raise argparse.ArgumentTypeError("must be greater than 0 and at most 100")
    return parsed


def _regex(value: str) -> str:
    try:
        re.compile(value)
    except re.error as exc:
        raise argparse.ArgumentTypeError(f"invalid regular expression: {exc}") from exc
    return value


def _add_output_options(parser: argparse.ArgumentParser, *, self_comparison: bool = False) -> None:
    parser.add_argument(
        "--output-dir", type=_output_path, default=Path.cwd(), help="Output directory"
    )
    parser.add_argument(
        "--prefix",
        type=_prefix,
        default="kffractbias",
        help="Output filename prefix; successful runs replace existing result files",
    )
    if self_comparison:
        parser.add_argument(
            "--name", default="self", help="Genome label used in metadata and plots"
        )
    else:
        parser.add_argument(
            "--target-name", default="target", help="Target label used in metadata and plots"
        )
        parser.add_argument(
            "--query-name", default="query", help="Query label used in metadata and plots"
        )
    parser.add_argument(
        "--window-size",
        type=_positive_int,
        default=100,
        help="Genes per complete sliding window; partial windows are omitted (default: 100)",
    )
    parser.add_argument(
        "--step-size", type=_positive_int, default=1, help="Genes advanced per window (default: 1)"
    )
    parser.add_argument(
        "--max-output-rows",
        type=_positive_int,
        help="Maximum combined gene and window TSV data rows; checked before writing tables",
    )
    parser.add_argument(
        "--denominator",
        choices=("all", "syntenic"),
        default="all",
        help="Use all target genes or only target genes with a syntenic match",
    )
    if self_comparison:
        parser.add_argument(
            "--seqids", action="append", default=[], help="Comma-separated genome sequences"
        )
    else:
        parser.add_argument(
            "--target-seqids", action="append", default=[], help="Comma-separated target sequences"
        )
        parser.add_argument(
            "--query-seqids",
            action="append",
            default=[],
            help="Comma-separated query sequences; includes unmatched selected sequences",
        )
    parser.add_argument(
        "--exclude-seqid-regex", type=_regex, default="", help="Regex for sequences to exclude"
    )
    parser.add_argument(
        "--include-unmatched-query-seqids",
        action="store_true",
        help="Include query sequences without retained synteny pairs in output tables",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Omit plots and remove existing plots for this prefix on success",
    )
    parser.add_argument(
        "--keep-failed-work",
        action="store_true",
        help="Retain this run's private staging directory and logs on failure",
    )


def _add_isoform_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--isoform-policy",
        choices=("error", "longest", "all"),
        default="error",
        help="Multiple CDS per gene locus: reject (default), select longest CDS, or count all identifiers",
    )


def _add_annotation_options(parser: argparse.ArgumentParser) -> None:
    _add_isoform_option(parser)
    parser.add_argument("--target-cds", required=True, type=_path, help="Target CDS FASTA")
    parser.add_argument("--target-gff", required=True, type=_path, help="Target GFF3/GTF")
    parser.add_argument("--query-cds", required=True, type=_path, help="Query CDS FASTA")
    parser.add_argument("--query-gff", required=True, type=_path, help="Query GFF3/GTF")
    parser.add_argument("--target-feature", help="Target GFF feature override, e.g. mRNA")
    parser.add_argument(
        "--target-attribute", help="Target GFF identifier attribute override, e.g. ID"
    )
    parser.add_argument("--query-feature", help="Query GFF feature override, e.g. mRNA")
    parser.add_argument(
        "--query-attribute", help="Query GFF identifier attribute override, e.g. ID"
    )
    parser.add_argument(
        "--minimum-mapping-fraction",
        type=_unit_interval,
        default=1.0,
        help="Minimum fraction of CDS identifiers that must map to each GFF (default: 1.0)",
    )


def _add_self_annotation_options(parser: argparse.ArgumentParser) -> None:
    _add_isoform_option(parser)
    parser.add_argument("--cds", required=True, type=_path, help="Genome CDS FASTA")
    parser.add_argument("--gff", required=True, type=_path, help="Genome GFF3/GTF")
    parser.add_argument("--feature", help="GFF feature override, e.g. mRNA")
    parser.add_argument("--attribute", help="GFF identifier attribute override, e.g. ID")
    parser.add_argument(
        "--minimum-mapping-fraction",
        type=_unit_interval,
        default=1.0,
        help="Minimum fraction of CDS identifiers that must map to the GFF (default: 1.0)",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="kffractbias",
        description="Offline fractionation-bias analysis",
    )
    parser.add_argument("--version", action="version", version=f"kfFractBias {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    calculate = subparsers.add_parser(
        "calculate",
        help="Calculate fractionation bias from precomputed synteny",
        description="Calculate fractionation bias from JCVI anchors or SynMap genomic-coordinate output.",
    )
    calculate.add_argument(
        "--synteny", required=True, type=_path, help="JCVI anchors or SynMap output"
    )
    calculate.add_argument(
        "--target-bed", required=True, type=_path, help="Target gene-coordinate BED"
    )
    calculate.add_argument(
        "--query-bed", required=True, type=_path, help="Query gene-coordinate BED"
    )
    calculate.add_argument("--format", choices=("auto", "jcvi", "synmap"), default="auto")
    calculate.add_argument(
        "--self",
        action="store_true",
        help="Reuse self-synteny anchors with identical BEDs; report within-genome retention",
    )
    _add_output_options(calculate)

    compare = subparsers.add_parser(
        "compare",
        help="Compare two annotations end-to-end",
        description="Generate local JCVI/QUOTA-ALIGN synteny and calculate fractionation bias.",
    )
    _add_annotation_options(compare)
    compare.add_argument(
        "--quota", required=True, type=validate_quota, help="Expected target:query depth, e.g. 1:2"
    )
    compare.add_argument(
        "--cpus", type=_positive_int, default=1, help="Threads for sequence alignment"
    )
    compare.add_argument("--cscore", type=_unit_interval, default=0.7, help="JCVI C-score cutoff")
    compare.add_argument("--aligner", choices=("last", "blast"), default="last")
    compare.add_argument(
        "--blast-task",
        choices=("blastn", "dc-megablast", "megablast"),
        help="BLAST+ search task; requires --aligner blast (default: blastn)",
    )
    compare.add_argument(
        "--force", action="store_true", help="Replace an existing synteny work directory"
    )
    _add_output_options(compare)

    selfcompare = subparsers.add_parser(
        "selfcompare",
        help="Measure within-genome self-synteny retention",
        description=(
            "Generate non-redundant JCVI self-synteny, apply symmetric QUOTA-ALIGN screening, "
            "and calculate within-genome retention profiles. This is not an outgroup-based fractionation analysis."
        ),
    )
    _add_self_annotation_options(selfcompare)
    selfcompare.add_argument(
        "--depth",
        required=True,
        type=_positive_int,
        help="Maximum expected homeologous block depth on each self-comparison axis",
    )
    selfcompare.add_argument(
        "--cpus", type=_positive_int, default=1, help="Threads for sequence alignment"
    )
    selfcompare.add_argument(
        "--cscore", type=_unit_interval, default=0.7, help="JCVI C-score cutoff"
    )
    selfcompare.add_argument("--aligner", choices=("last", "blast"), default="last")
    selfcompare.add_argument(
        "--blast-task",
        choices=("blastn", "dc-megablast", "megablast"),
        help="BLAST+ search task; requires --aligner blast (default: blastn)",
    )
    selfcompare.add_argument(
        "--self-hit-percent",
        type=_percent,
        default=98.0,
        help="JCVI percent-identity cutoff used to remove near-self hits (default: 98)",
    )
    selfcompare.add_argument(
        "--diagonal-bound",
        type=_positive_int,
        default=300,
        help="Minimum gene-rank distance between intrachromosomal self anchors (default: 300)",
    )
    selfcompare.add_argument(
        "--force", action="store_true", help="Replace an existing synteny work directory"
    )
    _add_output_options(selfcompare, self_comparison=True)

    validate = subparsers.add_parser(
        "validate",
        help="Validate CDS-to-GFF identifier mapping",
        description="Validate that both CDS FASTA files can be mapped to their annotations without running synteny.",
    )
    _add_annotation_options(validate)
    validate.add_argument(
        "--pairwise",
        action="store_true",
        help="Also require disjoint selected target/query IDs, as compare does",
    )

    subparsers.add_parser("formats", help="Describe supported synteny formats")
    subparsers.add_parser("version", help="Print the kfFractBias version")
    return parser


def _analysis_config(
    args: argparse.Namespace,
    *,
    synteny: Path,
    target_bed: Path,
    query_bed: Path,
    metadata: dict[str, Any] | None = None,
    additional_inputs: dict[str, Path] | None = None,
) -> AnalysisConfig:
    self_comparison = args.command == "selfcompare"
    return AnalysisConfig(
        synteny_path=synteny,
        synteny_format=args.format if hasattr(args, "format") else "jcvi",
        target_bed=target_bed,
        query_bed=query_bed,
        output_dir=args.output_dir,
        prefix=args.prefix,
        target_name=args.name if self_comparison else args.target_name,
        query_name=args.name if self_comparison else args.query_name,
        window_size=args.window_size,
        step_size=args.step_size,
        denominator=args.denominator,
        analysis_mode="self_synteny_retention"
        if self_comparison or getattr(args, "self", False)
        else "pairwise_fractionation_bias",
        target_seqids=tuple(args.seqids if self_comparison else args.target_seqids),
        query_seqids=tuple(args.seqids if self_comparison else args.query_seqids),
        exclude_seqid_regex=args.exclude_seqid_regex,
        include_unmatched_query_seqids=args.include_unmatched_query_seqids,
        make_plot=not args.no_plot,
        collect_rows=False,
        keep_failed_work=args.keep_failed_work,
        max_output_rows=args.max_output_rows,
        metadata=metadata or {},
        additional_inputs=additional_inputs or {},
    )


def _print_result(result: AnalysisResult) -> None:
    print(f"genes\t{result.genes_path}")
    print(f"windows\t{result.windows_path}")
    print(f"summary\t{result.summary_path}")
    if result.pdf_path:
        print(f"plot_pdf\t{result.pdf_path}")
    if result.png_path:
        print(f"plot_png\t{result.png_path}")


def command_calculate(args: argparse.Namespace) -> int:
    result = calculate_fractionation_bias(
        _analysis_config(
            args, synteny=args.synteny, target_bed=args.target_bed, query_bed=args.query_bed
        )
    )
    _print_result(result)
    return 0


def _check_feature_attribute_pairs(args: argparse.Namespace) -> None:
    for label in ("target", "query"):
        feature = getattr(args, f"{label}_feature")
        attribute = getattr(args, f"{label}_attribute")
        if (feature is None) != (attribute is None):
            raise ValueError(
                f"--{label}-feature and --{label}-attribute must be specified together"
            )


def _before_alignment(
    args: argparse.Namespace, run: RunContext
) -> Callable[[PreparedGenome, PreparedGenome], None]:
    def check(target: PreparedGenome, query: PreparedGenome) -> None:
        with run.stage("preflight"):
            run.capture_inputs(
                {
                    "target_bed": target.bed_path,
                    "query_bed": query.bed_path,
                    "prepared_target_cds": target.cds_path,
                    "prepared_query_cds": query.cds_path,
                }
            )
            run.verify_inputs()
            config = _analysis_config(
                args,
                synteny=run.work_dir / "pending.anchors",
                target_bed=target.bed_path,
                query_bed=query.bed_path,
            )
            estimates = preflight_analysis(config, target.mapping.genes, query.mapping.genes)
            if config.make_plot:
                from .plotting import preflight_plot

                preflight_plot()
            (run.work_dir / "preflight.json").write_text(
                json.dumps(estimates, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            print(
                f"kffractbias: dense output upper bound: {estimates['gene_table_rows_upper_bound']:,} gene rows, {estimates['window_table_rows_upper_bound']:,} window rows",
                file=sys.stderr,
            )

    return check


def command_compare(args: argparse.Namespace) -> int:
    _check_feature_attribute_pairs(args)
    if args.blast_task is not None and args.aligner != "blast":
        raise ValueError("--blast-task requires --aligner blast")
    source_inputs = {
        "source_target_cds": args.target_cds,
        "source_target_gff": args.target_gff,
        "source_query_cds": args.query_cds,
        "source_query_gff": args.query_gff,
    }
    with analysis_run(
        args.output_dir,
        args.prefix,
        source_inputs,
        synteny=True,
        force=args.force,
        keep_failed_work=args.keep_failed_work,
    ) as run:
        frozen = run.capture_inputs(source_inputs)
        with run.stage("synteny"):
            synteny = run_pairwise_synteny(
                target_cds=frozen["source_target_cds"],
                target_gff=frozen["source_target_gff"],
                query_cds=frozen["source_query_cds"],
                query_gff=frozen["source_query_gff"],
                work_dir=run.work_dir,
                quota=args.quota,
                cpus=args.cpus,
                cscore=args.cscore,
                aligner=args.aligner,
                target_feature=args.target_feature,
                target_attribute=args.target_attribute,
                query_feature=args.query_feature,
                query_attribute=args.query_attribute,
                minimum_mapping_fraction=args.minimum_mapping_fraction,
                isoform_policy=args.isoform_policy,
                before_alignment=_before_alignment(args, run),
                blast_task=args.blast_task or "blastn",
            )
        metadata = {
            "synteny_generation": {
                "method": "JCVI MCscan with QUOTA-ALIGN",
                "quota": synteny.quota,
                "command": list(synteny.command),
                "commands": [list(command) for command in synteny.commands],
                "blast_task": synteny.blast_task,
                "tool_versions": synteny.tool_versions,
                "target_gff_mapping": synteny.target.mapping.metadata(),
                "query_gff_mapping": synteny.query.mapping.metadata(),
            }
        }
        result = calculate_fractionation_bias(
            _analysis_config(
                args,
                synteny=synteny.anchors_path,
                target_bed=synteny.target.bed_path,
                query_bed=synteny.query.bed_path,
                metadata=metadata,
                additional_inputs={
                    "source_target_cds": args.target_cds,
                    "source_target_gff": args.target_gff,
                    "source_query_cds": args.query_cds,
                    "source_query_gff": args.query_gff,
                },
            ),
            run=run,
        )
    _print_result(result)
    return 0


def command_selfcompare(args: argparse.Namespace) -> int:
    if (args.feature is None) != (args.attribute is None):
        raise ValueError("--feature and --attribute must be specified together")
    if args.blast_task is not None and args.aligner != "blast":
        raise ValueError("--blast-task requires --aligner blast")
    source_inputs = {"source_cds": args.cds, "source_gff": args.gff}
    with analysis_run(
        args.output_dir,
        args.prefix,
        source_inputs,
        synteny=True,
        force=args.force,
        keep_failed_work=args.keep_failed_work,
    ) as run:
        frozen = run.capture_inputs(source_inputs)
        with run.stage("synteny"):
            synteny = run_self_synteny(
                cds=frozen["source_cds"],
                gff=frozen["source_gff"],
                work_dir=run.work_dir,
                depth=args.depth,
                cpus=args.cpus,
                cscore=args.cscore,
                aligner=args.aligner,
                feature=args.feature,
                attribute=args.attribute,
                minimum_mapping_fraction=args.minimum_mapping_fraction,
                isoform_policy=args.isoform_policy,
                self_hit_percent=args.self_hit_percent,
                intrachromosomal_diagonal_bound=args.diagonal_bound,
                before_alignment=_before_alignment(args, run),
                blast_task=args.blast_task or "blastn",
            )
        metadata = {
            "interpretation": (
                "Within-genome self-synteny retention asymmetry conditional on extant annotated genes; "
                "not an outgroup-based fractionation-bias estimate."
            ),
            "synteny_generation": {
                "method": "JCVI chromosome-aware self-synteny with shared-genome QUOTA-ALIGN",
                "depth": synteny.depth,
                "quota": synteny.quota,
                "commands": [list(command) for command in synteny.commands],
                "blast_task": synteny.blast_task,
                "tool_versions": synteny.tool_versions,
                "identity_and_mirror_handling": {
                    "jcvi_native_self_mode": True,
                    "chromosome_aware_scan": True,
                    "shared_genome_quota_constraints": True,
                    "self_hit_percent": synteny.self_hit_percent,
                    "intrachromosomal_diagonal_bound_genes": synteny.intrachromosomal_diagonal_bound,
                    "symmetric_quota_screen": True,
                },
                "gff_mapping": synteny.genome.mapping.metadata(),
            },
        }
        result = calculate_fractionation_bias(
            _analysis_config(
                args,
                synteny=synteny.anchors_path,
                target_bed=synteny.genome.bed_path,
                query_bed=synteny.genome.bed_path,
                metadata=metadata,
                additional_inputs={
                    "source_cds": args.cds,
                    "source_gff": args.gff,
                },
            ),
            run=run,
        )
    _print_result(result)
    return 0


def command_validate(args: argparse.Namespace) -> int:
    _check_feature_attribute_pairs(args)
    with tempfile.TemporaryDirectory(prefix="kffractbias_validate_") as temporary:
        target = prepare_genome(
            "target",
            args.target_cds,
            args.target_gff,
            temporary,
            feature=args.target_feature,
            attribute=args.target_attribute,
            minimum_mapping_fraction=args.minimum_mapping_fraction,
            isoform_policy=args.isoform_policy,
        )
        query = prepare_genome(
            "query",
            args.query_cds,
            args.query_gff,
            temporary,
            feature=args.query_feature,
            attribute=args.query_attribute,
            minimum_mapping_fraction=args.minimum_mapping_fraction,
            isoform_policy=args.isoform_policy,
        )
        if args.pairwise:
            validate_disjoint_identifiers(
                (gene.gene_id for gene in target.mapping.genes),
                (gene.gene_id for gene in query.mapping.genes),
                "CDS/GFF",
            )
    print(
        json.dumps(
            {
                "target": target.mapping.metadata(),
                "query": query.mapping.metadata(),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


def command_formats() -> int:
    print("jcvi\tJCVI .anchors files; either target-query orientation is accepted")
    print("synmap\tSynMap/DAGCHAINER genomic-coordinate output with || subfields")
    print("auto\tDetect SynMap by || delimiters; otherwise treat input as JCVI anchors")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    started = time.perf_counter()
    try:
        analysis_commands = {
            "calculate": command_calculate,
            "compare": command_compare,
            "selfcompare": command_selfcompare,
        }
        if args.command in analysis_commands:
            status = analysis_commands[args.command](args)
            if status == 0:
                print(
                    f"kffractbias: completed in {time.perf_counter() - started:.3f} seconds "
                    "(including commit and cleanup)",
                    file=sys.stderr,
                )
            return status
        if args.command == "validate":
            return command_validate(args)
        if args.command == "formats":
            return command_formats()
        if args.command == "version":
            print(f"kfFractBias {__version__}")
            return 0
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"kffractbias: error: {exc}", file=sys.stderr)
        return 2
    parser.error(f"unknown command: {args.command}")
    return 2

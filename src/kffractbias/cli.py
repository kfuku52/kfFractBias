from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .analysis import AnalysisConfig, calculate_fractionation_bias
from .jcvi import prepare_genome, run_pairwise_synteny, run_self_synteny, validate_quota


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
    parser.add_argument("--output-dir", type=_output_path, default=Path.cwd(), help="Output directory")
    parser.add_argument("--prefix", type=_prefix, default="kffractbias", help="Output filename prefix")
    if self_comparison:
        parser.add_argument("--name", default="self", help="Genome label used in metadata and plots")
    else:
        parser.add_argument("--target-name", default="target", help="Target label used in metadata and plots")
        parser.add_argument("--query-name", default="query", help="Query label used in metadata and plots")
    parser.add_argument("--window-size", type=_positive_int, default=100, help="Genes per sliding window (default: 100)")
    parser.add_argument("--step-size", type=_positive_int, default=1, help="Genes advanced per window (default: 1)")
    parser.add_argument(
        "--denominator",
        choices=("all", "syntenic"),
        default="all",
        help="Use all target genes or only target genes with a syntenic match",
    )
    if self_comparison:
        parser.add_argument("--seqids", action="append", default=[], help="Comma-separated genome sequences")
    else:
        parser.add_argument("--target-seqids", action="append", default=[], help="Comma-separated target sequences")
        parser.add_argument("--query-seqids", action="append", default=[], help="Comma-separated query sequences")
    parser.add_argument("--exclude-seqid-regex", type=_regex, default="", help="Regex for sequences to exclude")
    parser.add_argument("--no-plot", action="store_true", help="Do not create PDF and PNG plots")


def _add_annotation_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--target-cds", required=True, type=_path, help="Target CDS FASTA")
    parser.add_argument("--target-gff", required=True, type=_path, help="Target GFF3/GTF")
    parser.add_argument("--query-cds", required=True, type=_path, help="Query CDS FASTA")
    parser.add_argument("--query-gff", required=True, type=_path, help="Query GFF3/GTF")
    parser.add_argument("--target-feature", help="Target GFF feature override, e.g. mRNA")
    parser.add_argument("--target-attribute", help="Target GFF identifier attribute override, e.g. ID")
    parser.add_argument("--query-feature", help="Query GFF feature override, e.g. mRNA")
    parser.add_argument("--query-attribute", help="Query GFF identifier attribute override, e.g. ID")
    parser.add_argument(
        "--minimum-mapping-fraction",
        type=_unit_interval,
        default=0.5,
        help="Minimum fraction of CDS identifiers that must map to each GFF (default: 0.5)",
    )


def _add_self_annotation_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--cds", required=True, type=_path, help="Genome CDS FASTA")
    parser.add_argument("--gff", required=True, type=_path, help="Genome GFF3/GTF")
    parser.add_argument("--feature", help="GFF feature override, e.g. mRNA")
    parser.add_argument("--attribute", help="GFF identifier attribute override, e.g. ID")
    parser.add_argument(
        "--minimum-mapping-fraction",
        type=_unit_interval,
        default=0.5,
        help="Minimum fraction of CDS identifiers that must map to the GFF (default: 0.5)",
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
    calculate.add_argument("--synteny", required=True, type=_path, help="JCVI anchors or SynMap output")
    calculate.add_argument("--target-bed", required=True, type=_path, help="Target gene-coordinate BED")
    calculate.add_argument("--query-bed", required=True, type=_path, help="Query gene-coordinate BED")
    calculate.add_argument("--format", choices=("auto", "jcvi", "synmap"), default="auto")
    _add_output_options(calculate)

    compare = subparsers.add_parser(
        "compare",
        help="Compare two annotations end-to-end",
        description="Generate local JCVI/QUOTA-ALIGN synteny and calculate fractionation bias.",
    )
    _add_annotation_options(compare)
    compare.add_argument("--quota", required=True, type=validate_quota, help="Expected target:query depth, e.g. 1:2")
    compare.add_argument("--cpus", type=_positive_int, default=1, help="Threads for sequence alignment")
    compare.add_argument("--cscore", type=_unit_interval, default=0.7, help="JCVI C-score cutoff")
    compare.add_argument("--aligner", choices=("last", "blast"), default="last")
    compare.add_argument("--force", action="store_true", help="Replace an existing synteny work directory")
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
    selfcompare.add_argument("--cpus", type=_positive_int, default=1, help="Threads for sequence alignment")
    selfcompare.add_argument("--cscore", type=_unit_interval, default=0.7, help="JCVI C-score cutoff")
    selfcompare.add_argument("--aligner", choices=("last", "blast"), default="last")
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
        help="Minimum gene-rank distance for intrachromosomal self-synteny blocks (default: 300)",
    )
    selfcompare.add_argument("--force", action="store_true", help="Replace an existing synteny work directory")
    _add_output_options(selfcompare, self_comparison=True)

    validate = subparsers.add_parser(
        "validate",
        help="Validate CDS-to-GFF identifier mapping",
        description="Validate that both CDS FASTA files can be mapped to their annotations without running synteny.",
    )
    _add_annotation_options(validate)

    subparsers.add_parser("formats", help="Describe supported synteny formats")
    subparsers.add_parser("version", help="Print the kfFractBias version")
    return parser


def _analysis_config(
    args: argparse.Namespace,
    *,
    synteny: Path,
    target_bed: Path,
    query_bed: Path,
    metadata=None,
    additional_inputs=None,
):
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
        analysis_mode="self_synteny_retention" if self_comparison else "pairwise_fractionation_bias",
        target_seqids=tuple(args.seqids if self_comparison else args.target_seqids),
        query_seqids=tuple(args.seqids if self_comparison else args.query_seqids),
        exclude_seqid_regex=args.exclude_seqid_regex,
        make_plot=not args.no_plot,
        metadata=metadata or {},
        additional_inputs=additional_inputs or {},
    )


def _print_result(result) -> None:
    print(f"genes\t{result.genes_path}")
    print(f"windows\t{result.windows_path}")
    print(f"summary\t{result.summary_path}")
    if result.pdf_path:
        print(f"plot_pdf\t{result.pdf_path}")
    if result.png_path:
        print(f"plot_png\t{result.png_path}")


def command_calculate(args: argparse.Namespace) -> int:
    result = calculate_fractionation_bias(
        _analysis_config(args, synteny=args.synteny, target_bed=args.target_bed, query_bed=args.query_bed)
    )
    _print_result(result)
    return 0


def _check_feature_attribute_pairs(args: argparse.Namespace) -> None:
    for label in ("target", "query"):
        feature = getattr(args, f"{label}_feature")
        attribute = getattr(args, f"{label}_attribute")
        if (feature is None) != (attribute is None):
            raise ValueError(f"--{label}-feature and --{label}-attribute must be specified together")


def _prepare_synteny_work_dir(args: argparse.Namespace) -> Path:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = args.output_dir / f"{args.prefix}.synteny"
    if work_dir.exists():
        if not args.force:
            raise ValueError(f"Synteny work directory already exists: {work_dir}; use --force to replace it")
        if work_dir == work_dir.parent or work_dir.name in {"", ".", ".."}:
            raise ValueError(f"Refusing to replace unsafe work directory: {work_dir}")
        shutil.rmtree(work_dir)
    return work_dir


def command_compare(args: argparse.Namespace) -> int:
    _check_feature_attribute_pairs(args)
    work_dir = _prepare_synteny_work_dir(args)

    synteny = run_pairwise_synteny(
        target_cds=args.target_cds,
        target_gff=args.target_gff,
        query_cds=args.query_cds,
        query_gff=args.query_gff,
        work_dir=work_dir,
        quota=args.quota,
        cpus=args.cpus,
        cscore=args.cscore,
        aligner=args.aligner,
        target_feature=args.target_feature,
        target_attribute=args.target_attribute,
        query_feature=args.query_feature,
        query_attribute=args.query_attribute,
        minimum_mapping_fraction=args.minimum_mapping_fraction,
    )
    metadata = {
        "synteny_generation": {
            "method": "JCVI MCscan with QUOTA-ALIGN",
            "quota": synteny.quota,
            "command": list(synteny.command),
            "target_gff_mapping": {
                "feature": synteny.target.mapping.feature,
                "attribute": synteny.target.mapping.attribute,
                "matched_gene_count": synteny.target.mapping.matched_gene_count,
                "fasta_gene_count": synteny.target.mapping.fasta_gene_count,
            },
            "query_gff_mapping": {
                "feature": synteny.query.mapping.feature,
                "attribute": synteny.query.mapping.attribute,
                "matched_gene_count": synteny.query.mapping.matched_gene_count,
                "fasta_gene_count": synteny.query.mapping.fasta_gene_count,
            },
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
        )
    )
    _print_result(result)
    return 0


def command_selfcompare(args: argparse.Namespace) -> int:
    if (args.feature is None) != (args.attribute is None):
        raise ValueError("--feature and --attribute must be specified together")
    work_dir = _prepare_synteny_work_dir(args)
    synteny = run_self_synteny(
        cds=args.cds,
        gff=args.gff,
        work_dir=work_dir,
        depth=args.depth,
        cpus=args.cpus,
        cscore=args.cscore,
        aligner=args.aligner,
        feature=args.feature,
        attribute=args.attribute,
        minimum_mapping_fraction=args.minimum_mapping_fraction,
        self_hit_percent=args.self_hit_percent,
        intrachromosomal_diagonal_bound=args.diagonal_bound,
    )
    metadata = {
        "interpretation": (
            "Within-genome self-synteny retention asymmetry conditional on extant annotated genes; "
            "not an outgroup-based fractionation-bias estimate."
        ),
        "synteny_generation": {
            "method": "JCVI native self-synteny with symmetric QUOTA-ALIGN",
            "depth": synteny.depth,
            "quota": synteny.quota,
            "commands": [list(command) for command in synteny.commands],
            "identity_and_mirror_handling": {
                "jcvi_native_self_mode": True,
                "self_hit_percent": synteny.self_hit_percent,
                "intrachromosomal_diagonal_bound_genes": synteny.intrachromosomal_diagonal_bound,
                "symmetric_quota_screen": True,
            },
            "gff_mapping": {
                "feature": synteny.genome.mapping.feature,
                "attribute": synteny.genome.mapping.attribute,
                "matched_gene_count": synteny.genome.mapping.matched_gene_count,
                "fasta_gene_count": synteny.genome.mapping.fasta_gene_count,
            },
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
        )
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
        )
        query = prepare_genome(
            "query",
            args.query_cds,
            args.query_gff,
            temporary,
            feature=args.query_feature,
            attribute=args.query_attribute,
            minimum_mapping_fraction=args.minimum_mapping_fraction,
        )
    print(
        json.dumps(
            {
                "target": {
                    "feature": target.mapping.feature,
                    "attribute": target.mapping.attribute,
                    "matched_gene_count": target.mapping.matched_gene_count,
                    "fasta_gene_count": target.mapping.fasta_gene_count,
                },
                "query": {
                    "feature": query.mapping.feature,
                    "attribute": query.mapping.attribute,
                    "matched_gene_count": query.mapping.matched_gene_count,
                    "fasta_gene_count": query.mapping.fasta_gene_count,
                },
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
    try:
        if args.command == "calculate":
            return command_calculate(args)
        if args.command == "compare":
            return command_compare(args)
        if args.command == "selfcompare":
            return command_selfcompare(args)
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

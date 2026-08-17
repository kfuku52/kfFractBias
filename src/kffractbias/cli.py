from __future__ import annotations

import argparse
import json
import shutil
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .analysis import AnalysisConfig, calculate_fractionation_bias
from .jcvi import prepare_genome, run_pairwise_synteny, validate_quota


def _path(value: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"file not found: {value}")
    return path


def _output_path(value: str) -> Path:
    return Path(value).expanduser().resolve()


def _prefix(value: str) -> str:
    if not value or Path(value).name != value or value in {".", ".."}:
        raise argparse.ArgumentTypeError("prefix must be a non-empty filename component")
    return value


def _add_output_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--output-dir", type=_output_path, default=Path.cwd(), help="Output directory")
    parser.add_argument("--prefix", type=_prefix, default="kffractbias", help="Output filename prefix")
    parser.add_argument("--target-name", default="target", help="Target label used in metadata and plots")
    parser.add_argument("--query-name", default="query", help="Query label used in metadata and plots")
    parser.add_argument("--window-size", type=int, default=100, help="Genes per sliding window (default: 100)")
    parser.add_argument("--step-size", type=int, default=1, help="Genes advanced per window (default: 1)")
    parser.add_argument(
        "--denominator",
        choices=("all", "syntenic"),
        default="all",
        help="Use all target genes or only target genes with a syntenic match",
    )
    parser.add_argument("--target-seqids", action="append", default=[], help="Comma-separated target sequences")
    parser.add_argument("--query-seqids", action="append", default=[], help="Comma-separated query sequences")
    parser.add_argument("--exclude-seqid-regex", default="", help="Regex for sequences to exclude")
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
        type=float,
        default=0.5,
        help="Minimum fraction of CDS identifiers that must map to each GFF (default: 0.5)",
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
    compare.add_argument("--cpus", type=int, default=1, help="Threads for sequence alignment")
    compare.add_argument("--cscore", type=float, default=0.7, help="JCVI C-score cutoff")
    compare.add_argument("--aligner", choices=("last", "blast"), default="last")
    compare.add_argument("--force", action="store_true", help="Replace an existing synteny work directory")
    _add_output_options(compare)

    validate = subparsers.add_parser(
        "validate",
        help="Validate CDS-to-GFF identifier mapping",
        description="Validate that both CDS FASTA files can be mapped to their annotations without running synteny.",
    )
    _add_annotation_options(validate)

    subparsers.add_parser("formats", help="Describe supported synteny formats")
    subparsers.add_parser("version", help="Print the kfFractBias version")
    return parser


def _analysis_config(args: argparse.Namespace, *, synteny: Path, target_bed: Path, query_bed: Path, metadata=None):
    return AnalysisConfig(
        synteny_path=synteny,
        synteny_format=args.format if hasattr(args, "format") else "jcvi",
        target_bed=target_bed,
        query_bed=query_bed,
        output_dir=args.output_dir,
        prefix=args.prefix,
        target_name=args.target_name,
        query_name=args.query_name,
        window_size=args.window_size,
        step_size=args.step_size,
        denominator=args.denominator,
        target_seqids=tuple(args.target_seqids),
        query_seqids=tuple(args.query_seqids),
        exclude_seqid_regex=args.exclude_seqid_regex,
        make_plot=not args.no_plot,
        metadata=metadata or {},
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


def command_compare(args: argparse.Namespace) -> int:
    _check_feature_attribute_pairs(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = args.output_dir / f"{args.prefix}.synteny"
    if work_dir.exists():
        if not args.force:
            raise ValueError(f"Synteny work directory already exists: {work_dir}; use --force to replace it")
        if work_dir == work_dir.parent or work_dir.name in {"", ".", ".."}:
            raise ValueError(f"Refusing to replace unsafe work directory: {work_dir}")
        shutil.rmtree(work_dir)

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

"""Generate synthetic BED, anchor, and copy-mapping inputs for the README plots."""

from __future__ import annotations

import argparse
import math
import random
from pathlib import Path

GENES_PER_CHROMOSOME = 320
CHROMOSOMES = ("chr1", "chr2")
COPY_LABELS = ("A", "B")
QUERIES = {
    "A": ((0.83, 0.48), (0.72, 0.43)),
    "B": ((0.67, 0.61), (0.90, 0.70)),
}


def bed_row(chromosome: str, rank: int, gene: str) -> str:
    start = (rank - 1) * 1000
    return f"{chromosome}\t{start}\t{start + 900}\t{gene}\n"


def write_query(directory: Path, name: str, rates: tuple[tuple[float, ...], ...]) -> None:
    rng = random.Random(name)
    bed, anchors = [], []
    mapping = ["gene_id\tgroup_id\tcopy\tsegment_id\n"]
    for chrom, chromosome_rates in zip(CHROMOSOMES, rates, strict=True):
        for index, (copy, rate) in enumerate(zip(COPY_LABELS, chromosome_rates, strict=True)):
            seqid = f"{chrom}{copy}"
            anchors.append("###\n")
            for rank in range(1, GENES_PER_CHROMOSOME + 1):
                gene = f"{name}_{seqid}_{rank}"
                bed.append(bed_row(seqid, rank, gene))
                mapping.append(f"{gene}\t{chrom}\t{copy}\t{seqid}\n")
                # Vary detection along each copy and anchor both ends of the block.
                probability = rate + 0.10 * math.sin(rank / 24 + index * math.pi)
                if rank in (1, GENES_PER_CHROMOSOME) or rng.random() < probability:
                    anchors.append(f"ref_{chrom}_{rank}\t{gene}\t100\n")
    for suffix, rows in (("bed", bed), ("anchors", anchors), ("mapping.tsv", mapping)):
        (directory / f"{name}.{suffix}").write_text("".join(rows), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    directory = args.output_dir
    try:
        directory.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        parser.error("output directory already exists; choose a new directory")

    (directory / "reference.bed").write_text(
        "".join(
            bed_row(chrom, rank, f"ref_{chrom}_{rank}")
            for chrom in CHROMOSOMES
            for rank in range(1, GENES_PER_CHROMOSOME + 1)
        ),
        encoding="utf-8",
    )
    manifest = [
        "analysis_id\tspecies\texpected_copies\ttarget_bed\tquery_bed\tsynteny_file\tmapping_file\n"
    ]
    for label, rates in QUERIES.items():
        name = f"query-{label.lower()}"
        write_query(directory, name, rates)
        manifest.append(
            f"{name}\tQuery {label}\t2\treference.bed\t{name}.bed\t"
            f"{name}.anchors\t{name}.mapping.tsv\n"
        )
    (directory / "comparisons.tsv").write_text("".join(manifest), encoding="utf-8")


if __name__ == "__main__":
    main()

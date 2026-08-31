"""Generate reproducible synthetic CDS/GFF inputs for the annotation tutorial."""

from __future__ import annotations

import argparse
import random
from pathlib import Path


def gene_sequence(index: int, copy_index: int) -> str:
    generator = random.Random(index + 1729)
    sequence = list("ATG" + "".join(generator.choice("ACGT") for _ in range(300)) + "TAA")
    if copy_index:
        for position in range(10 + copy_index, len(sequence) - 3, 20):
            alternatives = "ACGT".replace(sequence[position], "")
            sequence[position] = alternatives[(index + position + copy_index) % len(alternatives)]
    return "".join(sequence)


def write_genome(directory: Path, name: str, seqids: tuple[str, ...]) -> None:
    fasta, annotation = [], ["##gff-version 3"]
    for copy_index, seqid in enumerate(seqids):
        for index in range(8):
            transcript = f"{name}{copy_index + 1}_{index + 1}"
            locus = f"{transcript}.gene"
            sequence = gene_sequence(index, copy_index)
            start = index * 1000 + 1
            end = start + len(sequence) - 1
            fasta.extend((f">{transcript}", sequence))
            for feature, attributes in (
                ("gene", f"ID={locus}"),
                ("mRNA", f"ID={transcript};Parent={locus}"),
                ("CDS", f"ID={transcript}.cds;Parent={transcript}"),
            ):
                phase = "0" if feature == "CDS" else "."
                annotation.append(
                    f"{seqid}\tsynthetic\t{feature}\t{start}\t{end}\t.\t+\t{phase}\t{attributes}"
                )
    (directory / f"{name}.cds.fa").write_text("\n".join(fasta) + "\n", encoding="utf-8")
    (directory / f"{name}.gff3").write_text("\n".join(annotation) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("annotation-inputs"))
    args = parser.parse_args()
    try:
        args.output_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        parser.error("output directory already exists; choose a new directory")
    write_genome(args.output_dir, "target", ("target_chr",))
    write_genome(args.output_dir, "query", ("query_a", "query_b"))
    print(f"Synthetic annotation inputs: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()

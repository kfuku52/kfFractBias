import hashlib
import json
import os
import random
from pathlib import Path

import pytest

from kffractbias.cli import main

pytestmark = pytest.mark.skipif(
    os.environ.get("KFFRACTBIAS_RUN_INTEGRATION") != "1",
    reason="set KFFRACTBIAS_RUN_INTEGRATION=1 to run JCVI/LAST integration",
)


def gene_sequence(index: int, copy_index: int = 0) -> str:
    generator = random.Random(index + 1729)
    sequence = list("ATG" + "".join(generator.choice("ACGT") for _ in range(300)) + "TAA")
    for position in range(10 + copy_index, len(sequence) - 3, 20):
        if copy_index:
            alternatives = "ACGT".replace(sequence[position], "")
            sequence[position] = alternatives[(index + position + copy_index) % len(alternatives)]
    return "".join(sequence)


def write_genome(tmp_path: Path, prefix: str, seqids: tuple[str, ...], copies: int):
    fasta_path = tmp_path / f"{prefix}.cds.fa"
    gff_path = tmp_path / f"{prefix}.gff3"
    fasta_lines = []
    gff_lines = ["##gff-version 3"]
    for copy_index in range(copies):
        seqid = seqids[copy_index]
        for gene_index in range(8):
            gene_id = f"{prefix}{copy_index + 1}_{gene_index + 1}"
            sequence = gene_sequence(gene_index, copy_index)
            fasta_lines.extend((f">{gene_id}", sequence))
            start = gene_index * 1000 + 1
            end = start + len(sequence) - 1
            gff_lines.append(f"{seqid}\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gene_id}")
    fasta_path.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")
    gff_path.write_text("\n".join(gff_lines) + "\n", encoding="utf-8")
    return fasta_path, gff_path


def test_compare_runs_jcvi_quota_align_offline(tmp_path):
    target_cds, target_gff = write_genome(tmp_path, "t", ("target_chr",), 1)
    query_cds, query_gff = write_genome(tmp_path, "q", ("query_a", "query_b"), 2)
    output_dir = tmp_path / "output"
    status = main(
        [
            "compare",
            "--target-cds",
            str(target_cds),
            "--target-gff",
            str(target_gff),
            "--query-cds",
            str(query_cds),
            "--query-gff",
            str(query_gff),
            "--quota",
            "1:2",
            "--window-size",
            "4",
            "--output-dir",
            str(output_dir),
            "--prefix",
            "synthetic",
        ]
    )
    assert status == 0
    assert (output_dir / "synthetic.synteny" / "target.query.lifted.1x2.anchors").is_file()
    assert (output_dir / "synthetic.genes.tsv").is_file()
    assert (output_dir / "synthetic.windows.tsv").is_file()
    assert (output_dir / "synthetic.summary.json").is_file()
    assert (output_dir / "synthetic.plot.pdf").is_file()
    assert (output_dir / "synthetic.plot.png").is_file()
    summary = json.loads((output_dir / "synthetic.summary.json").read_text(encoding="utf-8"))
    source_inputs = {
        "source_target_cds": target_cds,
        "source_target_gff": target_gff,
        "source_query_cds": query_cds,
        "source_query_gff": query_gff,
    }
    for label, path in source_inputs.items():
        assert summary["inputs"][label] == {
            "path": str(path.resolve()),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }


def test_selfcompare_runs_jcvi_quota_align_offline(tmp_path):
    cds, gff = write_genome(tmp_path, "s", ("self_a", "self_b"), 2)
    output_dir = tmp_path / "self-output"
    status = main(
        [
            "selfcompare",
            "--cds",
            str(cds),
            "--gff",
            str(gff),
            "--depth",
            "1",
            "--window-size",
            "4",
            "--output-dir",
            str(output_dir),
            "--prefix",
            "synthetic-self",
            "--no-plot",
        ]
    )
    assert status == 0
    assert (output_dir / "synthetic-self.synteny" / "self.self.lifted.1x1.anchors").is_file()
    summary = json.loads((output_dir / "synthetic-self.summary.json").read_text(encoding="utf-8"))
    assert summary["analysis_mode"] == "self_synteny_retention"
    assert summary["counts"]["synteny_pair_count"] > 0
    assert summary["metadata"]["synteny_generation"]["tool_versions"]["jcvi"]

import csv
import hashlib
import json
import os
import random
import shutil
from pathlib import Path

import pytest
from documentation import ANNOTATION_TUTORIAL, example_arguments, generate_annotation_inputs

from kffractbias.cli import build_parser, main
from kffractbias.io import parse_synteny_pairs, read_bed

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


@pytest.fixture(params=["last", "blast"])
def aligner(request):
    if request.param == "blast" and shutil.which("blastn") is None:
        pytest.skip("install BLAST+ to test the optional BLAST aligner")
    return request.param


def test_compare_runs_jcvi_quota_align_offline(tmp_path, aligner, monkeypatch):
    generate_annotation_inputs(tmp_path)
    monkeypatch.chdir(tmp_path)
    arguments = example_arguments(ANNOTATION_TUTORIAL, "compare")
    # Also exercise the tutorial's documented variant that enables plots.
    arguments.remove("--no-plot")
    arguments.extend(("--aligner", aligner, "--max-output-rows", "26"))
    options = build_parser().parse_args(arguments)
    target_cds, target_gff = options.target_cds, options.target_gff
    query_cds, query_gff = options.query_cds, options.query_gff
    output_dir, prefix = options.output_dir, options.prefix
    status = main(arguments)
    assert status == 0
    assert (output_dir / f"{prefix}.synteny" / "target.query.lifted.1x2.anchors").is_file()
    assert (output_dir / f"{prefix}.genes.tsv").is_file()
    assert (output_dir / f"{prefix}.windows.tsv").is_file()
    assert (output_dir / f"{prefix}.summary.json").is_file()
    assert (output_dir / f"{prefix}.plot.pdf").is_file()
    assert (output_dir / f"{prefix}.plot.png").is_file()
    summary = json.loads((output_dir / f"{prefix}.summary.json").read_text(encoding="utf-8"))
    expected = {(f"target1_{i}", f"query{copy}_{i}") for i in range(1, 9) for copy in (1, 2)}
    assert (
        set(
            parse_synteny_pairs(
                output_dir / f"{prefix}.synteny" / "target.query.lifted.1x2.anchors",
                "jcvi",
                {left for left, _ in expected},
                {right for _, right in expected},
            ).pairs
        )
        == expected
    )
    assert summary["counts"]["synteny_pair_count"] == 16
    assert summary["counts"]["gene_table_row_count"] == 16
    with (output_dir / f"{prefix}.windows.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 10
    assert all(row["retained_count"] == "4" and row["retention_fraction"] == "1" for row in rows)
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
    for label in (
        "target_bed",
        "query_bed",
        "prepared_target_cds",
        "prepared_query_cds",
        "synteny",
    ):
        path = Path(summary["inputs"][label]["path"])
        assert path.is_file()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == summary["inputs"][label]["sha256"]
    assert (output_dir / f"{prefix}.synteny" / "preflight.json").is_file()
    assert summary["metadata"]["synteny_generation"]["blast_task"] == (
        "blastn" if aligner == "blast" else None
    )


def test_selfcompare_runs_jcvi_quota_align_offline(tmp_path, aligner, monkeypatch):
    generate_annotation_inputs(tmp_path)
    monkeypatch.chdir(tmp_path)
    arguments = example_arguments(ANNOTATION_TUTORIAL, "selfcompare")
    arguments.extend(("--aligner", aligner, "--max-output-rows", "52"))
    options = build_parser().parse_args(arguments)
    output_dir, prefix = options.output_dir, options.prefix
    status = main(arguments)
    assert status == 0
    assert (output_dir / f"{prefix}.synteny" / "self.self.lifted.1x1.anchors").is_file()
    summary = json.loads((output_dir / f"{prefix}.summary.json").read_text(encoding="utf-8"))
    assert summary["analysis_mode"] == "self_synteny_retention"
    assert summary["counts"]["synteny_pair_count"] == 8
    assert summary["counts"]["interchromosomal_pair_count"] == 8
    assert summary["counts"]["intrachromosomal_pair_count"] == 0
    assert summary["counts"]["directed_synteny_pair_count"] == 16
    assert summary["counts"]["gene_table_row_count"] == 32
    assert summary["metadata"]["synteny_generation"]["tool_versions"]["jcvi"]
    expected = {(f"query1_{i}", f"query2_{i}") for i in range(1, 9)}
    identifiers = {gene for pair in expected for gene in pair}
    assert (
        set(
            parse_synteny_pairs(
                output_dir / f"{prefix}.synteny" / "self.self.lifted.1x1.anchors",
                "jcvi",
                identifiers,
                identifiers,
                allow_ambiguous_orientation=True,
            ).pairs
        )
        == expected
    )
    with (output_dir / f"{prefix}.windows.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 20
    assert all(
        row["retention_percent"] == ("0" if row["target_seqid"] == row["query_seqid"] else "100")
        for row in rows
    )


@pytest.mark.parametrize("depth,expected_count", [(1, 8), (2, 24)])
def test_real_self_quota_on_three_chromosomes(tmp_path, depth, expected_count):
    cds, gff = write_genome(tmp_path, "s", ("A", "B", "C"), 3)
    output = tmp_path / "out"
    assert (
        main(
            [
                "selfcompare",
                "--cds",
                str(cds),
                "--gff",
                str(gff),
                "--depth",
                str(depth),
                "--window-size",
                "4",
                "--output-dir",
                str(output),
                "--no-plot",
            ]
        )
        == 0
    )
    work = output / "kffractbias.synteny"
    identifiers = {gene.gene_id for gene in read_bed(work / "self.bed")}
    pairs = parse_synteny_pairs(
        work / f"self.self.lifted.{depth}x{depth}.anchors",
        "jcvi",
        identifiers,
        identifiers,
        allow_ambiguous_orientation=True,
    ).pairs
    assert len(pairs) == expected_count
    assert all(sum(gene in pair for pair in pairs) <= depth for gene in identifiers)

import csv
import json
import os
import subprocess
import sys
from dataclasses import replace

import pytest

from kffractbias.analysis import AnalysisConfig, calculate_fractionation_bias


def test_streamed_and_collected_tables_are_identical(tmp_path):
    target = tmp_path / "target.bed"
    query = tmp_path / "query.bed"
    anchors = tmp_path / "pairs.anchors"
    target.write_text("chr1\t0\t3\tt1\nchr1\t10\t13\tt2\nchr2\t0\t3\tt3\n")
    query.write_text("chrA\t0\t3\tq1\nchrB\t0\t3\tq2\n")
    anchors.write_text("t1\tq1\t10\nt3\tq2\t10\n")
    config = AnalysisConfig(
        anchors, "jcvi", target, query, tmp_path / "collected", window_size=2, make_plot=False
    )
    collected = calculate_fractionation_bias(config)
    streamed = calculate_fractionation_bias(
        replace(config, output_dir=tmp_path / "streamed", collect_rows=False)
    )
    assert streamed.gene_rows == streamed.window_rows == ()
    assert collected.genes_path.read_bytes() == streamed.genes_path.read_bytes()
    assert collected.windows_path.read_bytes() == streamed.windows_path.read_bytes()
    assert json.loads(streamed.summary_path.read_text())["counts"]["gene_table_row_count"] == len(
        collected.gene_rows
    )


def test_hash_seed_does_not_change_table_order(tmp_path):
    target = tmp_path / "target.bed"
    query = tmp_path / "query.bed"
    anchors = tmp_path / "pairs.anchors"
    target.write_text("chr1\t0\t3\tt1\n")
    query.write_text(
        "Q1\t0\t3\ta\nq1\t0\t3\tb\nq01\t0\t3\tc\nq1\t10\t13\tb01\nq1\t20\t23\tB1\nq1\t30\t33\tb1\n"
    )
    anchors.write_text("".join(f"t1\t{gene}\t10\n" for gene in ("a", "b", "c", "b01", "B1", "b1")))
    results = []
    for seed in (1, 2):
        output = tmp_path / str(seed)
        subprocess.run(
            [
                sys.executable,
                "-m",
                "kffractbias",
                "calculate",
                "--synteny",
                str(anchors),
                "--target-bed",
                str(target),
                "--query-bed",
                str(query),
                "--output-dir",
                str(output),
                "--window-size",
                "1",
                "--no-plot",
            ],
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
            check=True,
            capture_output=True,
        )
        results.append(
            (
                (output / "kffractbias.genes.tsv").read_bytes(),
                (output / "kffractbias.windows.tsv").read_bytes(),
            )
        )
    assert len(set(results)) == 1


@pytest.mark.parametrize(
    "denominator,expected",
    [
        (
            "all",
            [
                ("chrA", "t1", "t3", 1),
                ("chrA", "t3", "t5", 1),
                ("chrB", "t1", "t3", 1),
                ("chrB", "t3", "t5", 1),
            ],
        ),
        ("syntenic", [("chrA", "t1", "t4", 2), ("chrB", "t1", "t4", 1)]),
    ],
)
def test_retention_windows_count_genes_and_apply_step_after_filtering(
    tmp_path, denominator, expected
):
    target, query, anchors = (
        tmp_path / name for name in ("target.bed", "query.bed", "pairs.anchors")
    )
    # Unsorted coordinates, two unmatched genes, and a chromosome shorter than a window.
    target.write_text(
        "chr1\t50\t53\tt6\nchr1\t0\t3\tt1\nchr1\t40\t43\tt5\n"
        "chr1\t20\t23\tt3\nchr2\t0\t3\tt7\nchr1\t30\t33\tt4\nchr1\t10\t13\tt2\n"
    )
    query.write_text("chrA\t10\t13\tqa2\nchrB\t0\t3\tqb\nchrA\t0\t3\tqa1\n")
    # Multiple query hits still count as one retained target; reversed/duplicate pairs agree.
    anchors.write_text(
        "t1\tqa1\t10\nqa2\tt1\t10\nt1\tqa1\t10\nt3\tqb\t10\nt4\tqa1\t10\nqb\tt6\t10\nt7\tqa1\t10\n"
    )
    result = calculate_fractionation_bias(
        AnalysisConfig(
            anchors,
            "auto",
            target,
            query,
            tmp_path / "out",
            window_size=3,
            step_size=2,
            denominator=denominator,
            make_plot=False,
            collect_rows=False,
        )
    )
    with result.windows_path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["target_seqid"] for row in rows} == {"chr1"}
    assert [
        (row["query_seqid"], row["start_gene"], row["end_gene"], int(row["retained_count"]))
        for row in rows
    ] == expected
    assert [int(row["start_rank"]) for row in rows] == (
        [1, 3, 1, 3] if denominator == "all" else [1, 1]
    )
    assert [row["retention_fraction"] for row in rows] == (
        ["0.3333333333"] * 4 if denominator == "all" else ["0.6666666667", "0.3333333333"]
    )

import csv
import json
import os
import random
import subprocess
import sys
from dataclasses import replace

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
    for seed in range(1, 6):
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


def test_retention_matches_naive_oracle_for_120_cases(tmp_path):
    rng = random.Random(260831)
    for case in range(120):
        lengths = [rng.randrange(1, 31) for _ in range(rng.randrange(1, 5))]
        target_by_chr = {f"T{c}": [f"t{c}_{i}" for i in range(n)] for c, n in enumerate(lengths)}
        target_ids = [gene for genes in target_by_chr.values() for gene in genes]
        query_contigs = [f"Q{i}" for i in range(rng.randrange(1, 6))]
        query_to_chr = {
            f"q{j}_{i}": chrom for j, chrom in enumerate(query_contigs) for i in range(3)
        }
        pairs = {(t, q) for t in target_ids for q in query_to_chr if rng.random() < 0.11}
        pairs.add((target_ids[0], next(iter(query_to_chr))))
        target_rows = [
            f"{chrom}\t{i * 10}\t{i * 10 + 3}\t{gene}\n"
            for chrom, genes in target_by_chr.items()
            for i, gene in enumerate(genes)
        ]
        query_rows = [
            f"{chrom}\t{i * 10}\t{i * 10 + 3}\t{gene}\n"
            for i, (gene, chrom) in enumerate(query_to_chr.items())
        ]
        rng.shuffle(target_rows)
        rng.shuffle(query_rows)
        target, query, anchors = (
            tmp_path / name for name in ("target.bed", "query.bed", "pairs.anchors")
        )
        target.write_text("".join(target_rows))
        query.write_text("".join(query_rows))
        records = [(q, t) if rng.random() < 0.5 else (t, q) for t, q in sorted(pairs)]
        anchors.write_text(
            "".join(f"{left}\t{right}\t10\n" for left, right in records + records[:3])
        )
        window, step = rng.randrange(1, 12), rng.randrange(1, 7)
        denominator = rng.choice(["all", "syntenic"])
        result = calculate_fractionation_bias(
            AnalysisConfig(
                anchors,
                "auto",
                target,
                query,
                tmp_path / "out",
                window_size=window,
                step_size=step,
                denominator=denominator,
                include_unmatched_query_seqids=True,
                make_plot=False,
                collect_rows=False,
            )
        )
        expected = {}
        for chrom, genes in target_by_chr.items():
            if denominator == "syntenic":
                genes = [gene for gene in genes if any(t == gene for t, _ in pairs)]
            for qchrom in query_contigs:
                for start in range(0, len(genes) - window + 1, step):
                    selected = genes[start : start + window]
                    retained = sum(
                        any(t == gene and query_to_chr[q] == qchrom for t, q in pairs)
                        for gene in selected
                    )
                    expected[(chrom, qchrom, start + 1)] = (
                        retained,
                        selected[0],
                        selected[-1],
                        f"{retained / window:.10g}",
                    )
        with result.windows_path.open() as handle:
            actual = {
                (row["target_seqid"], row["query_seqid"], int(row["start_rank"])): (
                    int(row["retained_count"]),
                    row["start_gene"],
                    row["end_gene"],
                    row["retention_fraction"],
                )
                for row in csv.DictReader(handle, delimiter="\t")
            }
        assert actual == expected, f"oracle mismatch in case {case}"

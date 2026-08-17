import csv
import json
from pathlib import Path

from kffractbias.analysis import AnalysisConfig, calculate_fractionation_bias


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def make_inputs(tmp_path):
    target = write(
        tmp_path / "target.bed",
        "chr1\t0\t10\tt1\t0\t+\n"
        "chr1\t20\t30\tt2\t0\t+\n"
        "chr1\t40\t50\tt3\t0\t+\n"
        "chr1\t60\t70\tt4\t0\t+\n",
    )
    query = write(
        tmp_path / "query.bed",
        "chrA\t0\t10\tq1\t0\t+\n"
        "chrA\t20\t30\tq2\t0\t+\n"
        "chrB\t0\t10\tq3\t0\t+\n",
    )
    anchors = write(
        tmp_path / "target.query.anchors",
        "###\nt1\tq1\t100\nt2\tq2\t100\n###\nt3\tq3\t100\n",
    )
    return target, query, anchors


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_calculate_writes_expected_retention_windows(tmp_path):
    target, query, anchors = make_inputs(tmp_path)
    result = calculate_fractionation_bias(
        AnalysisConfig(
            synteny_path=anchors,
            synteny_format="jcvi",
            target_bed=target,
            query_bed=query,
            output_dir=tmp_path / "out",
            prefix="demo",
            window_size=2,
            make_plot=False,
        )
    )
    rows = read_tsv(result.windows_path)
    assert len(rows) == 6
    chr_a = [float(row["retention_percent"]) for row in rows if row["query_seqid"] == "chrA"]
    chr_b = [float(row["retention_percent"]) for row in rows if row["query_seqid"] == "chrB"]
    assert chr_a == [100.0, 50.0, 0.0]
    assert chr_b == [0.0, 50.0, 50.0]
    summary = json.loads(result.summary_path.read_text(encoding="utf-8"))
    assert summary["counts"]["synteny_pair_count"] == 3
    assert summary["parameters"]["denominator"] == "all"


def test_syntenic_denominator_removes_unmatched_target_genes(tmp_path):
    target, query, anchors = make_inputs(tmp_path)
    result = calculate_fractionation_bias(
        AnalysisConfig(
            synteny_path=anchors,
            synteny_format="auto",
            target_bed=target,
            query_bed=query,
            output_dir=tmp_path / "out",
            denominator="syntenic",
            window_size=2,
            make_plot=False,
        )
    )
    summary = json.loads(result.summary_path.read_text(encoding="utf-8"))
    assert summary["counts"]["target_gene_count"] == 4
    assert summary["counts"]["analyzed_target_gene_count"] == 3
    assert len(read_tsv(result.windows_path)) == 4


def test_exact_sequence_filtering(tmp_path):
    target, query, anchors = make_inputs(tmp_path)
    result = calculate_fractionation_bias(
        AnalysisConfig(
            synteny_path=anchors,
            synteny_format="jcvi",
            target_bed=target,
            query_bed=query,
            output_dir=tmp_path / "out",
            window_size=2,
            query_seqids=("chrA",),
            make_plot=False,
        )
    )
    assert {row["query_seqid"] for row in result.window_rows} == {"chrA"}


import csv
import hashlib
import json
from pathlib import Path

import pytest

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
    assert summary["parameters"]["synteny_format"] == "jcvi"
    assert summary["parameters"]["requested_synteny_format"] == "jcvi"


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
    assert summary["parameters"]["synteny_format"] == "jcvi"
    assert summary["parameters"]["requested_synteny_format"] == "auto"


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


def test_invalid_exclusion_regex_is_reported_as_value_error(tmp_path):
    target, query, anchors = make_inputs(tmp_path)
    with pytest.raises(ValueError, match="Invalid sequence exclusion regex"):
        calculate_fractionation_bias(
            AnalysisConfig(
                synteny_path=anchors,
                synteny_format="jcvi",
                target_bed=target,
                query_bed=query,
                output_dir=tmp_path / "out",
                exclude_seqid_regex="[",
                make_plot=False,
            )
        )


def test_additional_inputs_are_hashed_in_summary(tmp_path):
    target, query, anchors = make_inputs(tmp_path)
    source = write(tmp_path / "source.fa", ">gene1\nATG\n")
    result = calculate_fractionation_bias(
        AnalysisConfig(
            synteny_path=anchors,
            synteny_format="jcvi",
            target_bed=target,
            query_bed=query,
            output_dir=tmp_path / "out",
            window_size=2,
            make_plot=False,
            additional_inputs={"source_target_cds": source},
        )
    )
    summary = json.loads(result.summary_path.read_text(encoding="utf-8"))
    assert summary["inputs"]["source_target_cds"] == {
        "path": str(source),
        "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    }


def test_self_synteny_removes_identity_and_mirrors_then_maps_both_directions(tmp_path):
    bed = write(
        tmp_path / "self.bed",
        "chr1\t0\t10\tg1\n"
        "chr1\t20\t30\tg2\n"
        "chr2\t0\t10\tg3\n"
        "chr2\t20\t30\tg4\n",
    )
    anchors = write(
        tmp_path / "self.self.anchors",
        "###\ng1\tg1\t100\ng1\tg3\t90\ng3\tg1\t90\n###\ng2\tg4\t80\n",
    )
    result = calculate_fractionation_bias(
        AnalysisConfig(
            synteny_path=anchors,
            synteny_format="jcvi",
            target_bed=bed,
            query_bed=bed,
            output_dir=tmp_path / "out",
            analysis_mode="self_synteny_retention",
            window_size=1,
            make_plot=False,
        )
    )
    retained_pairs = {
        (row["target_gene"], row["query_genes"])
        for row in result.gene_rows
        if row["retained"] == 1
    }
    assert retained_pairs == {("g1", "g3"), ("g2", "g4"), ("g3", "g1"), ("g4", "g2")}
    summary = json.loads(result.summary_path.read_text(encoding="utf-8"))
    assert summary["analysis_mode"] == "self_synteny_retention"
    assert summary["counts"]["input_synteny_pair_count"] == 4
    assert summary["counts"]["synteny_pair_count"] == 2
    assert summary["counts"]["directed_synteny_pair_count"] == 4
    assert summary["counts"]["removed_identity_pair_count"] == 1
    assert summary["counts"]["removed_mirrored_pair_count"] == 1
    assert summary["counts"]["interchromosomal_pair_count"] == 2


def test_self_synteny_requires_identical_bed_coordinates(tmp_path):
    target = write(tmp_path / "target.bed", "chr1\t0\t10\tg1\nchr2\t0\t10\tg2\n")
    query = write(tmp_path / "query.bed", "chr1\t1\t10\tg1\nchr2\t0\t10\tg2\n")
    anchors = write(tmp_path / "self.anchors", "g1\tg2\t10\n")
    with pytest.raises(ValueError, match="identical target and query BED"):
        calculate_fractionation_bias(
            AnalysisConfig(
                synteny_path=anchors,
                synteny_format="jcvi",
                target_bed=target,
                query_bed=query,
                output_dir=tmp_path / "out",
                analysis_mode="self_synteny_retention",
                make_plot=False,
            )
        )

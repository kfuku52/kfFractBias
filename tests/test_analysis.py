import csv
import json
from dataclasses import replace
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
        "chrA\t0\t10\tq1\t0\t+\nchrA\t20\t30\tq2\t0\t+\nchrB\t0\t10\tq3\t0\t+\n",
    )
    anchors = write(
        tmp_path / "target.query.anchors",
        "###\nt1\tq1\t100\nt2\tq2\t100\n###\nt3\tq3\t100\n",
    )
    return target, query, anchors


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def test_self_synteny_removes_identity_and_mirrors_then_maps_both_directions(tmp_path):
    bed = write(
        tmp_path / "self.bed",
        "chr1\t0\t10\tg1\nchr1\t20\t30\tg2\nchr2\t0\t10\tg3\nchr2\t20\t30\tg4\n",
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
        (row["target_gene"], row["query_genes"]) for row in result.gene_rows if row["retained"] == 1
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


def test_output_limit_checks_exact_combined_rows_before_replacing_results(tmp_path, capsys):
    target, query, anchors = make_inputs(tmp_path)
    config = AnalysisConfig(
        anchors,
        "jcvi",
        target,
        query,
        tmp_path / "out",
        window_size=2,
        make_plot=False,
        max_output_rows=14,
    )
    result = calculate_fractionation_bias(config)
    before = {
        path: path.read_bytes()
        for path in (result.genes_path, result.windows_path, result.summary_path)
    }
    with pytest.raises(ValueError, match="requires 14 rows.*max-output-rows 13"):
        calculate_fractionation_bias(replace(config, max_output_rows=13))
    assert {path: path.read_bytes() for path in before} == before
    assert "8 gene rows, 6 window rows" in capsys.readouterr().err
    result = calculate_fractionation_bias(
        replace(config, denominator="syntenic", max_output_rows=10)
    )
    assert len(result.gene_rows) + len(result.window_rows) == 10
    assert json.loads(result.summary_path.read_text())["parameters"]["max_output_rows"] == 10

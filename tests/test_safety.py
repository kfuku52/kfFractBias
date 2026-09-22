import hashlib
import json
from pathlib import Path

import pytest

from kffractbias.analysis import AnalysisConfig, calculate_fractionation_bias
from kffractbias.cli import main


def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    target = write(tmp_path / "target.bed", "chr1\t0\t10\tt1\nchr1\t20\t30\tt2\n")
    query = write(tmp_path / "query.bed", "chrA\t0\t10\tq1\nchrB\t0\t10\tq2\n")
    anchors = write(tmp_path / "pairs.anchors", "t1\tq1\t10\n")
    return target, query, anchors


def test_input_output_collision_is_rejected_without_overwrite(tmp_path: Path) -> None:
    target, query, _anchors = inputs(tmp_path)
    output_dir = tmp_path / "output"
    synteny = write(output_dir / "demo.genes.tsv", "t1\tq1\t10\n")
    original = synteny.read_bytes()

    with pytest.raises(ValueError, match="same path"):
        calculate_fractionation_bias(
            AnalysisConfig(
                synteny_path=synteny,
                synteny_format="jcvi",
                target_bed=target,
                query_bed=query,
                output_dir=output_dir,
                prefix="demo",
                make_plot=False,
            )
        )

    assert synteny.read_bytes() == original


def test_pairwise_analysis_rejects_overlapping_gene_identifiers(tmp_path: Path) -> None:
    target = write(tmp_path / "target.bed", "chr1\t0\t10\tg1\n")
    query = write(tmp_path / "query.bed", "chrA\t0\t10\tg1\n")
    anchors = write(tmp_path / "pairs.anchors", "g1\tg1\t10\n")

    with pytest.raises(ValueError, match="must be disjoint"):
        calculate_fractionation_bias(
            AnalysisConfig(
                synteny_path=anchors,
                synteny_format="jcvi",
                target_bed=target,
                query_bed=query,
                output_dir=tmp_path / "output",
                make_plot=False,
            )
        )


def test_summary_records_hashes_runtime_and_parser_counts(tmp_path: Path) -> None:
    target, query, anchors = inputs(tmp_path)
    anchors.write_text("t1\tq1\t10\nt1\tq1\t10\n", encoding="utf-8")
    source = write(tmp_path / "source.fa", ">gene1\nATG\n")
    result = calculate_fractionation_bias(
        AnalysisConfig(
            synteny_path=anchors,
            synteny_format="jcvi",
            target_bed=target,
            query_bed=query,
            output_dir=tmp_path / "output",
            make_plot=False,
            additional_inputs={"source_target_cds": source},
        )
    )
    summary = json.loads(result.summary_path.read_text(encoding="utf-8"))

    assert summary["inputs"]["source_target_cds"] == {
        "path": str(source),
        "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    }
    assert summary["schema_version"] == 3
    assert summary["counts"]["input_synteny_record_count"] == 2
    assert summary["counts"]["duplicate_synteny_pair_count"] == 1
    assert summary["runtime"]["python"]
    source = summary["runtime"]["source"]
    assert source["hash_format"] == "kffractbias-python-source-v1"
    assert len(source["python_source_sha256"]) == 64
    assert (
        summary["output_sha256"]["genes"]
        == hashlib.sha256(result.genes_path.read_bytes()).hexdigest()
    )
    assert (
        summary["output_sha256"]["windows"]
        == hashlib.sha256(result.windows_path.read_bytes()).hexdigest()
    )


def test_unmatched_query_sequences_are_opt_in(tmp_path: Path) -> None:
    target, query, anchors = inputs(tmp_path)
    compact = calculate_fractionation_bias(
        AnalysisConfig(
            synteny_path=anchors,
            synteny_format="jcvi",
            target_bed=target,
            query_bed=query,
            output_dir=tmp_path / "compact",
            window_size=1,
            make_plot=False,
        )
    )
    complete = calculate_fractionation_bias(
        AnalysisConfig(
            synteny_path=anchors,
            synteny_format="jcvi",
            target_bed=target,
            query_bed=query,
            output_dir=tmp_path / "complete",
            window_size=1,
            include_unmatched_query_seqids=True,
            make_plot=False,
        )
    )

    assert {row["query_seqid"] for row in compact.gene_rows} == {"chrA"}
    assert {row["query_seqid"] for row in complete.gene_rows} == {"chrA", "chrB"}


def test_failed_replacement_restores_previous_synteny_directory(
    tmp_path: Path, monkeypatch
) -> None:
    target_cds = write(tmp_path / "target.fa", ">t1\nATG\n")
    target_gff = write(tmp_path / "target.gff", "chr1\tt\tmRNA\t1\t3\t.\t+\t.\tID=t1\n")
    query_cds = write(tmp_path / "query.fa", ">q1\nATG\n")
    query_gff = write(tmp_path / "query.gff", "chrA\tt\tmRNA\t1\t3\t.\t+\t.\tID=q1\n")
    output_dir = tmp_path / "output"
    sentinel = write(output_dir / "demo.synteny" / "sentinel.txt", "keep\n")

    def fail_synteny(**kwargs):
        work_dir = kwargs["work_dir"]
        work_dir.mkdir()
        write(work_dir / "partial.txt", "partial\n")
        raise RuntimeError("synthetic failure")

    monkeypatch.setattr("kffractbias.cli.run_pairwise_synteny", fail_synteny)
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
            "1:1",
            "--output-dir",
            str(output_dir),
            "--prefix",
            "demo",
            "--force",
        ]
    )

    assert status == 2
    assert sentinel.read_text(encoding="utf-8") == "keep\n"
    assert not (output_dir / "demo.synteny" / "partial.txt").exists()


def test_force_rejects_inputs_inside_replaceable_work_directory(
    tmp_path: Path, monkeypatch
) -> None:
    output_dir = tmp_path / "output"
    work_dir = output_dir / "demo.synteny"
    target_cds = write(work_dir / "target.fa", ">t1\nATG\n")
    target_gff = write(work_dir / "target.gff", "chr1\tt\tmRNA\t1\t3\t.\t+\t.\tID=t1\n")
    query_cds = write(tmp_path / "query.fa", ">q1\nATG\n")
    query_gff = write(tmp_path / "query.gff", "chrA\tt\tmRNA\t1\t3\t.\t+\t.\tID=q1\n")

    def unexpected_call(**_kwargs):
        raise AssertionError("synteny command must not run")

    monkeypatch.setattr("kffractbias.cli.run_pairwise_synteny", unexpected_call)
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
            "1:1",
            "--output-dir",
            str(output_dir),
            "--prefix",
            "demo",
            "--force",
        ]
    )

    assert status == 2
    assert target_cds.is_file()
    assert target_gff.is_file()

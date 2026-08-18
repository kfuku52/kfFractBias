from pathlib import Path

import pytest

from kffractbias.cli import build_parser, main
from kffractbias.io import AnnotationMapping, Gene
from kffractbias.jcvi import PreparedGenome, SelfSyntenyRun


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_public_subcommand_names_have_no_hyphens():
    parser = build_parser()
    subparsers = next(action for action in parser._actions if action.dest == "command")
    assert set(subparsers.choices) == {"calculate", "compare", "selfcompare", "validate", "formats", "version"}
    assert all("-" not in command for command in subparsers.choices)


def test_calculate_cli(tmp_path):
    target = write(tmp_path / "target.bed", "chr1\t0\t3\tt1\n")
    query = write(tmp_path / "query.bed", "chrA\t0\t3\tq1\n")
    anchors = write(tmp_path / "pairs.anchors", "t1\tq1\t10\n")
    output = tmp_path / "output"
    status = main(
        [
            "calculate",
            "--target-bed",
            str(target),
            "--query-bed",
            str(query),
            "--synteny",
            str(anchors),
            "--window-size",
            "1",
            "--output-dir",
            str(output),
            "--no-plot",
        ]
    )
    assert status == 0
    assert (output / "kffractbias.genes.tsv").is_file()
    assert (output / "kffractbias.windows.tsv").is_file()
    assert (output / "kffractbias.summary.json").is_file()


def test_selfcompare_cli_uses_one_annotation_and_self_mode(tmp_path, monkeypatch):
    cds = write(tmp_path / "genome.cds.fa", ">g1\nATG\n>g2\nATG\n")
    gff = write(
        tmp_path / "genome.gff3",
        "##gff-version 3\nchr1\ttest\tmRNA\t1\t3\t.\t+\t.\tID=g1\n"
        "chr2\ttest\tmRNA\t1\t3\t.\t+\t.\tID=g2\n",
    )

    def fake_self_synteny(**kwargs):
        work_dir = kwargs["work_dir"]
        work_dir.mkdir(parents=True)
        bed = write(work_dir / "self.bed", "chr1\t0\t3\tg1\nchr2\t0\t3\tg2\n")
        anchors = write(work_dir / "self.self.lifted.1x1.anchors", "###\ng1\tg2\t10\n")
        mapping = AnnotationMapping(
            feature="mRNA",
            attribute="ID",
            genes=(Gene("chr1", 0, 3, "g1", "+"), Gene("chr2", 0, 3, "g2", "+")),
            fasta_gene_count=2,
            matched_gene_count=2,
        )
        genome = PreparedGenome("self", bed, work_dir / "self.cds", mapping)
        return SelfSyntenyRun(
            genome=genome,
            anchors_path=anchors,
            commands=(("python", "ortholog"), ("python", "quota")),
            depth=1,
            quota="1:1",
            self_hit_percent=98.0,
            intrachromosomal_diagonal_bound=300,
        )

    monkeypatch.setattr("kffractbias.cli.run_self_synteny", fake_self_synteny)
    output = tmp_path / "output"
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
            "1",
            "--output-dir",
            str(output),
            "--no-plot",
        ]
    )
    assert status == 0
    summary = (output / "kffractbias.summary.json").read_text(encoding="utf-8")
    assert '"analysis_mode": "self_synteny_retention"' in summary
    assert '"synteny_pair_count": 1' in summary


def test_formats_cli(capsys):
    assert main(["formats"]) == 0
    output = capsys.readouterr().out
    assert "jcvi" in output
    assert "synmap" in output


@pytest.mark.parametrize(
    ("option", "value"),
    (
        ("--window-size", "0"),
        ("--step-size", "-1"),
        ("--exclude-seqid-regex", "["),
    ),
)
def test_calculate_rejects_invalid_analysis_options_during_parsing(option, value):
    parser = build_parser()
    with pytest.raises(SystemExit) as caught:
        parser.parse_args(
            [
                "calculate",
                "--synteny",
                __file__,
                "--target-bed",
                __file__,
                "--query-bed",
                __file__,
                option,
                value,
            ]
        )
    assert caught.value.code == 2

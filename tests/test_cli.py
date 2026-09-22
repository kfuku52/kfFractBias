from pathlib import Path

import pytest

from kffractbias.cli import build_parser, main


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_calculate_reuses_self_anchors_without_alignment(tmp_path):
    import json

    bed = write(tmp_path / "self.bed", "chr1\t0\t3\tg1\nchr2\t0\t3\tg2\n")
    anchors = write(tmp_path / "self.anchors", "g1\tg2\t10\ng2\tg1\t10\n")
    output = tmp_path / "out"
    assert (
        main(
            [
                "calculate",
                "--self",
                "--target-bed",
                str(bed),
                "--query-bed",
                str(bed),
                "--synteny",
                str(anchors),
                "--output-dir",
                str(output),
                "--window-size",
                "1",
                "--no-plot",
            ]
        )
        == 0
    )
    summary = json.loads((output / "kffractbias.summary.json").read_text())
    assert summary["analysis_mode"] == "self_synteny_retention"
    assert summary["counts"]["synteny_pair_count"] == 1
    assert summary["counts"]["directed_synteny_pair_count"] == 2
    assert "not an outgroup-based" in summary["metadata"]["interpretation"]


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

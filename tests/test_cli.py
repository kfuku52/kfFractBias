from pathlib import Path

from kffractbias.cli import build_parser, main


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_public_subcommand_names_have_no_hyphens():
    parser = build_parser()
    subparsers = next(action for action in parser._actions if action.dest == "command")
    assert set(subparsers.choices) == {"calculate", "compare", "validate", "formats", "version"}
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


def test_formats_cli(capsys):
    assert main(["formats"]) == 0
    output = capsys.readouterr().out
    assert "jcvi" in output
    assert "synmap" in output


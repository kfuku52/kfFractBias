import json
import re
import shlex
import sys
from pathlib import Path

import pytest
from test_run import compare_args, compare_inputs

from kffractbias.cli import build_parser, main
from kffractbias.jcvi import _run_checked, preflight_tools


def test_readme_commands_parse_without_external_tools(monkeypatch):
    monkeypatch.setattr("kffractbias.cli._path", Path)
    parser = build_parser()
    readme = (Path(__file__).resolve().parents[1] / "README.md").read_text()
    commands = re.findall(r"^kffractbias .+$", readme.replace("\\\n", " "), flags=re.MULTILINE)
    assert len(commands) >= 5
    for command in commands:
        parser.parse_args(shlex.split(command.replace("\\\n", " "))[1:])


def test_public_self_defaults_and_options():
    args = build_parser().parse_args(
        ["selfcompare", "--cds", __file__, "--gff", __file__, "--depth", "1"]
    )
    assert (args.diagonal_bound, args.self_hit_percent, args.isoform_policy) == (300, 98.0, "error")
    assert not args.keep_failed_work


def test_blast_task_override_requires_blast(tmp_path, capsys):
    args = compare_args(compare_inputs(tmp_path), tmp_path / "out")
    assert main(args + ["--blast-task", "megablast"]) == 2
    assert "--blast-task requires --aligner blast" in capsys.readouterr().err


@pytest.mark.parametrize("extra", [["--target-seqids", "missing"], ["--query-seqids", "missing"]])
def test_bad_sequence_selection_fails_before_alignment(tmp_path, monkeypatch, extra, capsys):
    def unexpected(*args, **kwargs):
        pytest.fail("external processing must not start")

    monkeypatch.setattr("kffractbias.jcvi._run_checked", unexpected)
    monkeypatch.setattr("kffractbias.jcvi.preflight_tools", unexpected)
    paths = compare_inputs(tmp_path)
    assert main(compare_args(paths, tmp_path / "out") + extra) == 2
    assert "sequence identifiers" in capsys.readouterr().err


def test_plot_dependency_is_checked_before_alignment(tmp_path, monkeypatch):
    def missing():
        raise RuntimeError("matplotlib unavailable")

    def unexpected(*args, **kwargs):
        pytest.fail("external processing must not start")

    monkeypatch.setattr("kffractbias.plotting.preflight_plot", missing)
    monkeypatch.setattr("kffractbias.jcvi._run_checked", unexpected)
    args = compare_args(compare_inputs(tmp_path), tmp_path / "out")
    args.remove("--no-plot")
    assert main(args) == 2


@pytest.mark.parametrize(
    "aligner,expected", [("last", "lastal, lastdb"), ("blast", "blastn, makeblastdb")]
)
def test_missing_aligner_reports_required_executables(monkeypatch, aligner, expected):
    monkeypatch.setattr("kffractbias.jcvi.shutil.which", lambda _name: None)
    with pytest.raises(RuntimeError, match=expected):
        preflight_tools(aligner)


def test_command_failure_records_stderr_command_and_timing(tmp_path):
    command = (
        sys.executable,
        "-c",
        "import sys; print('diagnostic', file=sys.stderr); sys.exit(7)",
    )
    with pytest.raises(RuntimeError, match="diagnostic"):
        _run_checked(command, cwd=tmp_path, description="failing stage")
    log = next((tmp_path / "logs").glob("*.log"))
    record = json.loads(next((tmp_path / "logs").glob("*.json")).read_text())
    assert "diagnostic" in log.read_text()
    assert record["command"] == list(command)
    assert record["returncode"] == 7
    assert record["elapsed_seconds"] > 0


def test_empty_fasta_header_is_a_cli_input_error(tmp_path, capsys):
    paths = compare_inputs(tmp_path)
    paths["target_cds"].write_text(">\nATG\n")
    args = ["validate"]
    for label, path in paths.items():
        args.extend(("--" + label.replace("_", "-"), str(path)))
    assert main(args) == 2
    assert "Empty FASTA identifier" in capsys.readouterr().err

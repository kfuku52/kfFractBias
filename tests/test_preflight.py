import json
import sys

import pytest
from test_run import compare_args, compare_inputs

from kffractbias.cli import main
from kffractbias.jcvi import _run_checked, preflight_tools


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


@pytest.mark.parametrize("command", ["validate", "compare", "selfcompare"])
def test_invalid_gff_is_rejected_before_alignment(tmp_path, monkeypatch, capsys, command):
    def unexpected(*args, **kwargs):
        pytest.fail("external alignment must not start for invalid annotation")

    monkeypatch.setattr("kffractbias.jcvi.preflight_tools", unexpected)
    monkeypatch.setattr("kffractbias.jcvi._run_checked", unexpected)
    paths = compare_inputs(tmp_path)
    paths["target_gff"].write_text("chr1\ttest\tmRNA\t1\t3\t.\tINVALID\t.\tID=t1\n")
    output = tmp_path / "out"
    if command == "compare":
        args = compare_args(paths, output)
    elif command == "selfcompare":
        args = [
            command,
            "--cds",
            str(paths["target_cds"]),
            "--gff",
            str(paths["target_gff"]),
            "--depth",
            "1",
            "--output-dir",
            str(output),
            "--no-plot",
        ]
    else:
        args = [command]
        for label, path in paths.items():
            args.extend(("--" + label.replace("_", "-"), str(path)))
    assert main(args) == 2
    assert "Invalid GFF strand" in capsys.readouterr().err
    assert not list(output.glob("*.summary.json"))


def test_pairwise_validation_rejects_overlap_but_mapping_only_still_allows_self(tmp_path, capsys):
    paths = compare_inputs(tmp_path)
    args = ["validate"]
    for axis in ("target", "query"):
        for kind in ("cds", "gff"):
            args.extend((f"--{axis}-{kind}", str(paths[f"target_{kind}"])))
    assert main(args) == 0
    assert main([*args, "--pairwise"]) == 2
    assert "must be disjoint" in capsys.readouterr().err
    args = ["validate", "--pairwise"]
    for label, path in paths.items():
        args.extend(("--" + label.replace("_", "-"), str(path)))
    assert main(args) == 0

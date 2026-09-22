import json
import multiprocessing
import os
from pathlib import Path

import pytest

from kffractbias.cli import main
from kffractbias.run import RecoveryRequired, analysis_run


def compare_inputs(tmp_path):
    paths = {}
    for label, gene, chrom in (("target", "t1", "chr1"), ("query", "q1", "chrA")):
        cds = tmp_path / f"{label}.fa"
        cds.write_text(f">{gene}\nATG\n")
        gff = tmp_path / f"{label}.gff"
        gff.write_text(f"{chrom}\ttest\tmRNA\t1\t3\t.\t+\t.\tID={gene}\n")
        paths[f"{label}_cds"], paths[f"{label}_gff"] = cds, gff
    return paths


def compare_args(paths, output):
    args = [
        "compare",
        "--quota",
        "1:1",
        "--no-plot",
        "--output-dir",
        str(output),
        "--prefix",
        "demo",
    ]
    for label, path in paths.items():
        args.extend(("--" + label.replace("_", "-"), str(path)))
    return args


def hold_run(output, ready, release):
    with analysis_run(Path(output), "demo", {}, synteny=True, force=True) as run:
        run.work_dir.mkdir()
        marker = run.work_dir / "running.txt"
        marker.write_text("owned by first process")
        ready.put(str(marker))
        if not release.wait(15):
            raise RuntimeError("test parent did not release worker")


@pytest.mark.parametrize("command", ["compare", "selfcompare", "calculate"])
def test_all_commands_hold_the_same_lock_before_touching_work(tmp_path, command, capsys):
    paths = compare_inputs(tmp_path)
    output = tmp_path / "output"
    old = output / "demo.synteny" / "old.txt"
    old.parent.mkdir(parents=True)
    old.write_text("old work")
    context = multiprocessing.get_context("spawn")
    ready, release = context.Queue(), context.Event()
    process = context.Process(target=hold_run, args=(str(output), ready, release))
    process.start()
    try:
        marker = Path(ready.get(timeout=15))
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
                "--prefix",
                "demo",
                "--no-plot",
            ]
        else:
            bed = tmp_path / "target.bed"
            bed.write_text("chr1\t0\t3\tt1\n")
            query = tmp_path / "query.bed"
            query.write_text("chrA\t0\t3\tq1\n")
            pairs = tmp_path / "pairs.anchors"
            pairs.write_text("t1\tq1\t10\n")
            args = [
                command,
                "--synteny",
                str(pairs),
                "--target-bed",
                str(bed),
                "--query-bed",
                str(query),
                "--output-dir",
                str(output),
                "--prefix",
                "demo",
                "--no-plot",
            ]
        if command != "calculate":
            args.append("--force")
        assert main(args) == 2
        assert "Another analysis" in capsys.readouterr().err
        assert marker.read_text() == "owned by first process"
        assert old.read_text() == "old work"
    finally:
        release.set()
        process.join(15)
        if process.is_alive():
            process.terminate()
            process.join(5)
    assert process.exitcode == 0


@pytest.mark.parametrize("link_type", ["symlink", "hardlink"])
def test_lock_links_do_not_modify_unrelated_files(tmp_path, link_type):
    sentinel = tmp_path / "sentinel"
    sentinel.write_text("do not overwrite")
    output = tmp_path / "output"
    output.mkdir()
    lock = output / ".demo.lock"
    if link_type == "symlink":
        lock.symlink_to(sentinel)
    else:
        os.link(sentinel, lock)
    with pytest.raises((OSError, ValueError)):
        with analysis_run(output, "demo", {}):
            pytest.fail("unsafe lock acquired")
    assert sentinel.read_text() == "do not overwrite"


def test_input_snapshot_rejects_changes_during_copy(tmp_path, monkeypatch):
    from kffractbias import run as run_module

    source = tmp_path / "source"
    source.write_text("before")
    copy = run_module.shutil.copyfile

    def mutate(source_path, destination):
        copy(source_path, destination)
        source.write_text("after")

    monkeypatch.setattr(run_module.shutil, "copyfile", mutate)
    with pytest.raises(RuntimeError, match="while taking a snapshot"):
        with analysis_run(tmp_path / "out", "demo", {"source": source}):
            pytest.fail("changing source accepted")
    assert not list((tmp_path / "out").glob(".demo.staging-*"))


def test_compare_uses_frozen_source_and_rejects_late_mutation(tmp_path, monkeypatch, capsys):
    paths = compare_inputs(tmp_path)
    output = tmp_path / "output"
    old = output / "demo.synteny" / "old.txt"
    old.parent.mkdir(parents=True)
    old.write_text("previous work")
    summary = output / "demo.summary.json"
    summary.write_text("previous summary")
    copied = []

    def fake_command(command, *, cwd, description):
        copied.append((cwd / "target.cds").read_text())
        paths["target_cds"].write_text(">t1\nCCC\n")
        (cwd / "target.query.lifted.1x1.anchors").write_text("t1\tq1\t10\n")

    monkeypatch.setattr("kffractbias.jcvi._run_checked", fake_command)
    monkeypatch.setattr("kffractbias.jcvi.preflight_tools", lambda _aligner: None)
    assert main(compare_args(paths, output) + ["--force", "--keep-failed-work"]) == 2
    assert copied == [">t1\nATG\n"]
    assert old.read_text() == "previous work"
    assert summary.read_text() == "previous summary"
    assert "Input files changed during analysis" in capsys.readouterr().err
    retained = list(output.glob(".demo.staging-*"))
    assert len(retained) == 1
    assert (retained[0] / "demo.synteny" / "target.cds").read_text() == copied[0]
    failure = json.loads((retained[0] / "failure.json").read_text())
    assert failure["inputs"]["source_target_cds"]["path"] == str(paths["target_cds"].resolve())


@pytest.mark.parametrize("operation", ["backup", "install"])
def test_commit_failure_restores_work_and_all_outputs(tmp_path, monkeypatch, operation):
    from kffractbias import run as run_module

    output = tmp_path / "out"
    output.mkdir()
    old_work = output / "demo.synteny"
    old_work.mkdir()
    (old_work / "marker").write_text("old work")
    for name in ("genes.tsv", "summary.json"):
        (output / f"demo.{name}").write_text(f"old {name}")
    replace = run_module.os.replace
    failed = False
    final_summary = output / "demo.summary.json"

    def fail_once(source, destination):
        nonlocal failed
        # Fail after work and the gene table have been installed, so rollback
        # must restore both completed outputs and the interrupted summary.
        affected = source if operation == "backup" else destination
        if not failed and affected == final_summary:
            failed = True
            raise OSError("injected commit failure")
        replace(source, destination)

    with pytest.raises(OSError, match="injected commit"):
        with analysis_run(output, "demo", {}, synteny=True, force=True) as run:
            run.work_dir.mkdir()
            (run.work_dir / "marker").write_text("new work")
            staged = {}
            for name in ("genes.tsv", "summary.json"):
                path = run.staging_dir / name
                path.write_text("new")
                staged[output / f"demo.{name}"] = path
            monkeypatch.setattr(run_module.os, "replace", fail_once)
            run.commit(staged)
    assert failed
    assert (old_work / "marker").read_text() == "old work"
    for name in ("genes.tsv", "summary.json"):
        assert (output / f"demo.{name}").read_text() == f"old {name}"
    assert not list(output.glob(".demo.staging-*"))


def test_failed_rollback_keeps_recovery_files(tmp_path, monkeypatch):
    from kffractbias import run as run_module

    output = tmp_path / "out"
    output.mkdir()
    final = output / "demo.genes.tsv"
    final.write_text("old")
    replace = run_module.os.replace

    def fail_install_and_restore(source, destination):
        if destination == final:
            raise OSError("filesystem unavailable")
        replace(source, destination)

    with pytest.raises(RecoveryRequired):
        with analysis_run(output, "demo", {}) as run:
            staged = run.staging_dir / "new.tsv"
            staged.write_text("new")
            monkeypatch.setattr(run_module.os, "replace", fail_install_and_restore)
            run.commit({final: staged})
    recovery = next(output.glob(".demo.staging-*"))
    assert next(recovery.glob("backup-*")).read_text() == "old"
    assert (recovery / "failure.json").is_file()


def test_force_cannot_remove_source_through_case_alias(tmp_path):
    output = tmp_path / "out"
    work = output / "demo.synteny"
    work.mkdir(parents=True)
    source = work / "source.fa"
    source.write_text(">x\nATG\n")
    if not (output / "DEMO.synteny").exists():
        pytest.skip("requires a case-insensitive filesystem")
    with pytest.raises(ValueError, match="inside the replaceable"):
        with analysis_run(output, "DEMO", {"source": source}, synteny=True, force=True):
            pytest.fail("would replace an input-containing directory")
    assert source.read_text() == ">x\nATG\n"

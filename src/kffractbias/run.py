"""Own one locked run, immutable input copies, and recoverable output commits."""

from __future__ import annotations

import fcntl
import json
import os
import shutil
import stat
import sys
import tempfile
import time
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .io import sha256_file


def validate_prefix(prefix: str) -> None:
    if (
        not prefix
        or Path(prefix).name != prefix
        or prefix in {".", ".."}
        or any(ord(character) < 32 or ord(character) == 127 for character in prefix)
    ):
        raise ValueError("prefix must be a non-empty filename component")


def output_paths(output_dir: Path, prefix: str) -> dict[str, Path]:
    validate_prefix(prefix)
    directory = output_dir.resolve()
    return {
        "genes": directory / f"{prefix}.genes.tsv",
        "windows": directory / f"{prefix}.windows.tsv",
        "summary": directory / f"{prefix}.summary.json",
        "plot_pdf": directory / f"{prefix}.plot.pdf",
        "plot_png": directory / f"{prefix}.plot.png",
        "lock": directory / f".{prefix}.lock",
    }


def validate_separation(inputs: Mapping[str, Path], outputs: Mapping[str, Path]) -> None:
    for input_label, input_path in inputs.items():
        for output_label, output_path in outputs.items():
            same = input_path.resolve() == output_path.resolve()
            try:
                same = same or input_path.samefile(output_path)
            except OSError:
                pass
            if same:
                raise ValueError(
                    f"Input {input_label!r} and output {output_label!r} refer to the same path: {input_path}"
                )


def _inside_directory(path: Path, directory: Path) -> bool:
    path, directory = path.resolve(), directory.resolve()
    if path.is_relative_to(directory):
        return True
    # resolve() does not normalize spelling on case-insensitive filesystems.
    for parent in path.parents:
        try:
            if parent.samefile(directory):
                return True
        except OSError:
            pass
    return False


@contextmanager
def exclusive_output_lock(lock_path: Path) -> Iterator[None]:
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    # Do not follow links, open special files in blocking mode, or truncate an
    # existing inode. The file remains in place after unlock to avoid split locks.
    descriptor = os.open(lock_path, os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW | os.O_NONBLOCK, 0o600)
    try:
        info = os.fstat(descriptor)
        if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
            raise ValueError(f"Lock must be a regular file without hard links: {lock_path}")
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(
                f"Another analysis is writing the same output prefix: {lock_path}"
            ) from exc
        current = lock_path.lstat()
        if (info.st_dev, info.st_ino) != (current.st_dev, current.st_ino):
            raise RuntimeError(f"Lock file was replaced while acquiring it: {lock_path}")
        yield
    finally:
        os.close(descriptor)


class RecoveryRequired(RuntimeError):
    """A rollback failed; never automatically delete the recovery directory."""


def commit_outputs(staged: dict[Path, Path | None], staging_dir: Path) -> None:
    backups: dict[Path, Path] = {}
    installed: dict[Path, Path] = {}
    for final_path, staged_path in staged.items():
        if final_path.is_symlink():
            raise ValueError(f"Refusing to replace symlinked output: {final_path}")
        if final_path.exists() and final_path.is_dir() != bool(
            staged_path and staged_path.is_dir()
        ):
            raise ValueError(f"Output path has the wrong file type: {final_path}")
    try:
        for index, (final_path, staged_path) in enumerate(staged.items()):
            if final_path.exists():
                backup = staging_dir / f"backup-{index}-{final_path.name}"
                os.replace(final_path, backup)
                backups[final_path] = backup
            if staged_path is not None:
                os.replace(staged_path, final_path)
                installed[final_path] = staged_path
    except BaseException as exc:
        errors: list[str] = []
        for final_path, staged_path in reversed(installed.items()):
            try:
                os.replace(final_path, staged_path)
            except OSError as rollback_error:
                errors.append(str(rollback_error))
        for final_path, backup in backups.items():
            try:
                os.replace(backup, final_path)
            except OSError as rollback_error:
                errors.append(str(rollback_error))
        if errors:
            raise RecoveryRequired(
                f"Output rollback failed; recovery files are in {staging_dir}: {'; '.join(errors)}"
            ) from exc
        raise


@dataclass
class RunContext:
    output_dir: Path
    prefix: str
    staging_dir: Path
    synteny: bool = False
    inputs: dict[str, dict[str, str]] = field(default_factory=dict)
    timings: dict[str, float] = field(default_factory=dict)
    _reads: dict[str, Path] = field(default_factory=dict)
    _checks: dict[Path, str] = field(default_factory=dict)
    _active: bool = field(default=False, init=False)
    _committed: bool = field(default=False, init=False)

    @property
    def work_dir(self) -> Path:
        return self.staging_dir / f"{self.prefix}.synteny"

    @property
    def final_work_dir(self) -> Path:
        return self.output_dir / self.work_dir.name

    def published_path(self, path: Path) -> Path:
        if self.synteny and path.is_relative_to(self.work_dir):
            return self.final_work_dir / path.relative_to(self.work_dir)
        return path

    def published_metadata(self, value: Any) -> Any:
        if isinstance(value, dict):
            return {key: self.published_metadata(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [self.published_metadata(item) for item in value]
        if isinstance(value, str) and self.synteny:
            return value.replace(str(self.work_dir), str(self.final_work_dir))
        return value

    def capture_inputs(self, inputs: Mapping[str, Path]) -> dict[str, Path]:
        """Freeze outside sources; watch already-private generated inputs."""
        self.require_active(self.output_dir, self.prefix)
        for label, source in inputs.items():
            source = source.resolve()
            if label in self.inputs:
                if self.inputs[label]["path"] != str(self.published_path(source)):
                    raise ValueError(f"Input label was reused for a different path: {label}")
                continue
            if not source.is_file():
                raise ValueError(f"Input must be a regular file: {source}")
            expected = sha256_file(source)
            if source.is_relative_to(self.staging_dir):
                frozen = source
            else:
                directory = self.staging_dir / "inputs"
                directory.mkdir(exist_ok=True)
                frozen = directory / f"{len(self.inputs)}-{source.name}"
                shutil.copyfile(source, frozen)
                if sha256_file(frozen) != expected or sha256_file(source) != expected:
                    raise RuntimeError(f"Input files changed while taking a snapshot: {label}")
            self.inputs[label] = {"path": str(self.published_path(source)), "sha256": expected}
            self._reads[label] = frozen
            self._checks[source] = expected
            self._checks[frozen] = expected
        return {label: self._reads[label] for label in inputs}

    def verify_inputs(self) -> None:
        for path, expected in self._checks.items():
            try:
                unchanged = sha256_file(path) == expected
            except OSError:
                unchanged = False
            if not unchanged:
                raise RuntimeError(f"Input files changed during analysis: {path}")

    @contextmanager
    def stage(self, name: str) -> Iterator[None]:
        start = time.perf_counter()
        try:
            yield
        finally:
            self.timings[name] = self.timings.get(name, 0.0) + time.perf_counter() - start

    def require_active(self, output_dir: Path, prefix: str) -> None:
        if (
            not self._active
            or self._committed
            or self.output_dir != output_dir.resolve()
            or self.prefix != prefix
        ):
            raise RuntimeError("Analysis must use its active, matching run context")

    def commit(self, staged: dict[Path, Path | None]) -> None:
        self.require_active(self.output_dir, self.prefix)
        self.verify_inputs()
        if self.synteny:
            staged = {self.final_work_dir: self.work_dir, **staged}
        with self.stage("commit"):
            commit_outputs(staged, self.staging_dir)
        self._committed = True


@contextmanager
def analysis_run(
    output_dir: Path,
    prefix: str,
    sources: Mapping[str, Path],
    *,
    synteny: bool = False,
    force: bool = False,
    keep_failed_work: bool = False,
) -> Iterator[RunContext]:
    outputs = output_paths(output_dir, prefix)
    output_dir = outputs["genes"].parent
    validate_separation(sources, outputs)
    with exclusive_output_lock(outputs["lock"]):
        for label, path in outputs.items():
            if label != "lock" and (path.is_symlink() or (path.exists() and not path.is_file())):
                raise ValueError(
                    f"Output must be a regular file, not a directory or symlink: {path}"
                )
        work = output_dir / f"{prefix}.synteny"
        if synteny:
            if work.is_symlink() or (work.exists() and not work.is_dir()):
                raise ValueError(f"Refusing to replace unsafe synteny work directory: {work}")
            if any(_inside_directory(path, work) for path in sources.values()):
                raise ValueError(
                    "Input files must not be located inside the replaceable synteny work directory"
                )
            if work.exists() and not force:
                raise ValueError(
                    f"Synteny work directory already exists: {work}; use --force to replace it"
                )
        staging = Path(tempfile.mkdtemp(prefix=f".{prefix}.staging-", dir=output_dir))
        run = RunContext(output_dir, prefix, staging, synteny)
        run._active = True
        retain = False
        try:
            with run.stage("snapshot"):
                run.capture_inputs(sources)
            yield run
        except BaseException as exc:
            retain = keep_failed_work or isinstance(exc, RecoveryRequired)
            if retain:
                # Keep this run's directory under its unique name, including any
                # rollback backups. Never rename or remove another run's work.
                try:
                    (staging / "failure.json").write_text(
                        json.dumps(
                            {
                                "error": str(exc),
                                "inputs": run.inputs,
                                "timings_seconds": run.timings,
                            },
                            indent=2,
                            sort_keys=True,
                        )
                        + "\n",
                        encoding="utf-8",
                    )
                except OSError as log_error:
                    print(f"kffractbias: could not write failure log: {log_error}", file=sys.stderr)
                print(f"kffractbias: failed work retained at {staging}", file=sys.stderr)
            raise
        finally:
            run._active = False
            if not retain:
                try:
                    shutil.rmtree(staging)
                except OSError as cleanup_error:
                    # A cleanup failure must not turn a successful commit into
                    # an apparent failed analysis, or mask its original error.
                    print(
                        f"kffractbias: could not remove staging directory {staging}: {cleanup_error}",
                        file=sys.stderr,
                    )

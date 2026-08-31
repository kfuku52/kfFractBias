"""Identify the package sources, independently of the caller's working directory."""

from __future__ import annotations

import hashlib
import os
import subprocess
from pathlib import Path
from typing import TypedDict

HASH_FORMAT = "kffractbias-python-source-v1"


class SourceMetadata(TypedDict):
    git_commit: str | None
    git_dirty: bool | None
    python_source_sha256: str | None
    hash_format: str


def _source_hash(package: Path) -> str | None:
    digest = hashlib.sha256(HASH_FORMAT.encode("ascii") + b"\0")
    try:
        sources = sorted(
            package.rglob("*.py"), key=lambda path: path.relative_to(package).as_posix()
        )
        if not sources:
            return None
        for path in sources:
            digest.update(path.relative_to(package).as_posix().encode("utf-8") + b"\0")
            digest.update(hashlib.sha256(path.read_bytes()).digest())
    except OSError:
        return None
    return digest.hexdigest()


def _git_state(package: Path) -> tuple[str | None, bool | None]:
    # An inherited GIT_DIR/GIT_WORK_TREE/index must not identify another checkout.
    environment = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
    environment["GIT_OPTIONAL_LOCKS"] = "0"

    def git(*args: str) -> str:
        return subprocess.check_output(  # nosec B603
            ("git", "-C", str(package), "--literal-pathspecs", *args),
            env=environment,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=2,
        ).strip()

    try:
        root = Path(git("rev-parse", "--show-toplevel")).resolve()
        relative = package.relative_to(root) / "__init__.py"
        # A wheel in an untracked/ignored venv is not the enclosing repository.
        git("ls-files", "--error-unmatch", "--", str(root / relative))
        commit = git("rev-parse", "--verify", "HEAD")
        dirty = bool(git("status", "--porcelain=v1", "--untracked-files=normal"))
        return commit, dirty
    except (OSError, ValueError, subprocess.SubprocessError):
        return None, None


def capture_source_metadata(package: Path | None = None) -> SourceMetadata:
    """Hash Python source bytes; Git is optional and never required to run."""
    package = (package or Path(__file__).parent).resolve()
    commit, dirty = _git_state(package)
    return {
        "git_commit": commit,
        "git_dirty": dirty,
        "python_source_sha256": _source_hash(package),
        "hash_format": HASH_FORMAT,
    }


# Keep the identity from import time, not a checkout edited during alignment.
# This describes files on disk at import, not dynamically patched Python objects.
_SOURCE = capture_source_metadata()


def source_metadata() -> SourceMetadata:
    return _SOURCE.copy()

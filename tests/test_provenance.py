import hashlib
import subprocess

import pytest

from kffractbias.provenance import HASH_FORMAT, capture_source_metadata


def git(root, *args):
    return subprocess.check_output(
        ["git", "-C", str(root), *args], text=True, stderr=subprocess.DEVNULL
    ).strip()


def checkout(root):
    package = root / "src" / "kffractbias"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text('"""A package fixture."""\n')
    (package / "calculation.py").write_text("value = 1\n")
    git(root, "init", "-q")
    git(root, "add", "src")
    git(
        root,
        "-c",
        "user.name=Test",
        "-c",
        "user.email=test@example.invalid",
        "-c",
        "commit.gpgsign=false",
        "commit",
        "-qm",
        f"fixture {root.name}",
    )
    return package


def test_source_identity_uses_the_package_checkout_not_cwd_or_git_environment(
    tmp_path, monkeypatch
):
    package = checkout(tmp_path / "package")
    foreign = tmp_path / "foreign"
    checkout(foreign)
    monkeypatch.chdir(foreign)
    monkeypatch.setenv("GIT_DIR", str(foreign / ".git"))
    monkeypatch.setenv("GIT_WORK_TREE", str(foreign))
    metadata = capture_source_metadata(package)
    expected = subprocess.check_output(
        ["git", "--git-dir", str(tmp_path / "package" / ".git"), "rev-parse", "HEAD"], text=True
    ).strip()
    assert metadata["git_commit"] == expected
    assert metadata["git_commit"] != git(foreign, "rev-parse", "HEAD")
    assert metadata["git_dirty"] is False
    (package / "calculation.py").write_text("value = 2\n")
    changed = capture_source_metadata(package)
    assert changed["git_commit"] == expected
    assert changed["git_dirty"] is True
    assert changed["python_source_sha256"] != metadata["python_source_sha256"]


def test_untracked_package_does_not_claim_enclosing_repository_commit(tmp_path):
    root = tmp_path / "unrelated"
    checkout(root)
    package = root / ".venv" / "lib" / "kffractbias"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("value = 1\n")
    metadata = capture_source_metadata(package)
    assert metadata["git_commit"] is None
    assert metadata["git_dirty"] is None
    assert metadata["python_source_sha256"]


def test_hash_covers_relative_names_and_source_bytes_but_not_bytecode(tmp_path):
    package = tmp_path / "package"
    package.mkdir()
    contents = {"__init__.py": b"value = 1\n", "sub/algorithm.py": b"value = 2\n"}
    expected = hashlib.sha256(HASH_FORMAT.encode("ascii") + b"\0")
    for name, content in contents.items():
        path = package / name
        path.parent.mkdir(exist_ok=True)
        path.write_bytes(content)
        expected.update(name.encode("utf-8") + b"\0" + hashlib.sha256(content).digest())
    original = capture_source_metadata(package)
    assert original["python_source_sha256"] == expected.hexdigest()
    cache = package / "__pycache__"
    cache.mkdir()
    (cache / "__init__.cpython-312.pyc").write_bytes(b"different bytecode")
    assert capture_source_metadata(package) == original
    (package / "sub" / "algorithm.py").rename(package / "sub" / "renamed.py")
    assert (
        capture_source_metadata(package)["python_source_sha256"] != original["python_source_sha256"]
    )


@pytest.mark.parametrize("failure", [FileNotFoundError(), subprocess.TimeoutExpired("git", 2)])
def test_git_unavailable_preserves_source_hash(tmp_path, monkeypatch, failure):
    from kffractbias import provenance

    (tmp_path / "__init__.py").write_text("value = 1\n")

    def unavailable(*args, **kwargs):
        raise failure

    monkeypatch.setattr(provenance.subprocess, "check_output", unavailable)
    metadata = capture_source_metadata(tmp_path)
    assert metadata["git_commit"] is None
    assert metadata["git_dirty"] is None
    assert len(metadata["python_source_sha256"]) == 64

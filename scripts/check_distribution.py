"""Test built wheel without dependencies and the sdist in a separate environment."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import tomllib
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    version = tomllib.loads((root / "pyproject.toml").read_text())["project"]["version"]
    (wheel,) = (root / "dist").glob(f"kffractbias-{version}-*.whl")
    archive = root / "dist" / f"kffractbias-{version}.tar.gz"
    environment = {
        key: value
        for key, value in os.environ.items()
        if key
        not in {
            "PYTHONPATH",
            "VIRTUAL_ENV",
            "UV_PROJECT_ENVIRONMENT",
            "KFFRACTBIAS_RUN_INTEGRATION",
        }
    }
    environment["PYTHONNOUSERSITE"] = "1"
    with tempfile.TemporaryDirectory(prefix="kffractbias-dist-") as temporary:
        work = Path(temporary)

        def run(command: list[str], cwd: Path = work) -> None:
            subprocess.run(command, cwd=cwd, env=environment, check=True)  # nosec B603

        venv = work / "wheel-env"
        run(["uv", "venv", "--python", sys.executable, str(venv)])
        python = str(venv / "bin" / "python")
        run(["uv", "pip", "install", "--python", python, "--no-deps", "--no-index", str(wheel)])
        run(
            [
                python,
                "-I",
                "-c",
                "import importlib.util, pathlib, sys, kffractbias; assert pathlib.Path(kffractbias.__file__).is_relative_to(sys.prefix); assert importlib.util.find_spec('matplotlib') is None",
            ]
        )
        example = work / "example"
        shutil.copytree(root / "examples" / "minimal", example)
        run(
            [
                python,
                "-I",
                "-m",
                "kffractbias",
                "calculate",
                "--target-bed",
                str(example / "target.bed"),
                "--query-bed",
                str(example / "query.bed"),
                "--synteny",
                str(example / "pairs.anchors"),
                "--window-size",
                "2",
                "--output-dir",
                str(work / "result"),
                "--no-plot",
            ]
        )
        summary = json.loads((work / "result" / "kffractbias.summary.json").read_text())
        assert summary["counts"]["synteny_pair_count"] > 0
        assert summary["program_version"] == version
        with tarfile.open(archive) as source:
            source.extractall(work / "source", filter="data")
        (extracted,) = (work / "source").iterdir()
        for path in (
            "CITATION.cff",
            "CONTRIBUTING.md",
            "CHANGELOG.md",
            "uv.lock",
            "scripts/check.py",
            "scripts/generate_example.py",
            "tests/test_metadata.py",
            "examples/minimal/target.bed",
            "examples/annotations/README.md",
            "docs/README.md",
            "docs/formats.md",
            "docs/methods.md",
            "docs/migration.md",
            "docs/troubleshooting.md",
        ):
            assert (extracted / path).is_file(), f"sdist missing {path}"
        run(
            [
                "uv",
                "sync",
                "--locked",
                "--extra",
                "test",
                "--extra",
                "plot",
                "--python",
                sys.executable,
            ],
            cwd=extracted,
        )
        run(
            [str(extracted / ".venv" / "bin" / "python"), "-I", "-m", "pytest", "-q"], cwd=extracted
        )
    print("Wheel and sdist passed isolated checks.")


if __name__ == "__main__":
    main()

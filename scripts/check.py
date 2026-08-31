"""Run fast local checks, or opt into package and real alignment checks."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full", action="store_true", help="Add coverage, build and isolated distribution tests"
    )
    parser.add_argument(
        "--integration", action="store_true", help="Require compare extra and LAST/BLAST on PATH"
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    commands = [
        ["-m", "ruff", "check", "src", "tests", "scripts"],
        ["-m", "ruff", "format", "--check", "src", "tests", "scripts"],
        ["-m", "mypy", "src"],
    ]
    commands += (
        [["-m", "coverage", "run", "-m", "pytest", "-q"]] if args.full else [["-m", "pytest", "-q"]]
    )
    if args.full:
        commands += [
            ["-m", "coverage", "report"],
            ["-m", "build"],
            ["scripts/check_distribution.py"],
        ]
    environment = os.environ.copy()
    if args.integration:
        environment["KFFRACTBIAS_RUN_INTEGRATION"] = "1"
    for command in commands:
        print("Running:", sys.executable, *command, flush=True)
        subprocess.run([sys.executable, *command], cwd=root, env=environment, check=True)  # nosec B603


if __name__ == "__main__":
    main()

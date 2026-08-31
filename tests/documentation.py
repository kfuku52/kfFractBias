"""Read runnable examples from the maintained Markdown, without running a shell."""

import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from markdown_it import MarkdownIt

ROOT = Path(__file__).resolve().parents[1]
ANNOTATION_TUTORIAL = ROOT / "examples" / "annotations" / "README.md"
MINIMAL_TUTORIAL = ROOT / "examples" / "minimal" / "README.md"
MARKDOWN = MarkdownIt("commonmark")
DOCUMENTS = sorted(
    [*ROOT.glob("*.md"), *(ROOT / "docs").rglob("*.md"), *(ROOT / "examples").rglob("*.md")]
)


@dataclass(frozen=True)
class ExampleCommand:
    line: int
    argv: tuple[str, ...]


def shell_commands(document):
    commands = []
    for token in MARKDOWN.parse(document.read_text(encoding="utf-8")):
        if token.type != "fence" or token.info.strip() not in {"bash", "sh", "shell", "console"}:
            continue
        pending = ""
        first_line = 0
        for offset, line in enumerate(token.content.splitlines()):
            if not pending:
                first_line = token.map[0] + offset + 2
            pending += line.strip().removeprefix("$ ")
            if pending.endswith("\\"):
                pending = pending[:-1] + " "
                continue
            argv = shlex.split(pending, comments=True)
            if argv:
                commands.append(ExampleCommand(first_line, tuple(argv)))
            pending = ""
        assert not pending, f"Unfinished shell continuation at {document}:{first_line}"
    return commands


def cli_commands(document):
    commands = []
    for command in shell_commands(document):
        words = command.argv
        if words[0] == "kffractbias":
            arguments = words[1:]
        elif words[:2] == ("uv", "run") and "kffractbias" in words:
            arguments = words[words.index("kffractbias") + 1 :]
        elif words[0] in {"python", "python3"} and words[1:3] == ("-m", "kffractbias"):
            arguments = words[3:]
        else:
            continue
        commands.append(ExampleCommand(command.line, arguments))
    return commands


def example_arguments(document, subcommand):
    (command,) = [command for command in cli_commands(document) if command.argv[0] == subcommand]
    return list(command.argv)


def generate_annotation_inputs(work):
    (command,) = [
        command
        for command in shell_commands(ANNOTATION_TUTORIAL)
        if command.argv[:2] == ("python", "scripts/generate_example.py")
    ]
    # Only this known synthetic generator is executed; documentation is not a shell script.
    subprocess.run(
        [sys.executable, str(ROOT / command.argv[1]), *command.argv[2:]],
        cwd=work,
        check=True,
        capture_output=True,
        text=True,
    )

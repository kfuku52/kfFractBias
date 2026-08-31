import re
import tomllib
from pathlib import Path

from kffractbias import __version__

ROOT = Path(__file__).resolve().parents[1]


def _declared_version(path: Path) -> str:
    match = re.search(
        r"^version\s*[:=]\s*[\"']?([^\"'\s]+)", path.read_text(encoding="utf-8"), re.MULTILINE
    )
    assert match is not None
    return match.group(1)


def test_version_metadata_stays_synchronized() -> None:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    expected = project["version"]
    assert _declared_version(ROOT / "CITATION.cff") == expected
    assert __version__ == expected
    lock = tomllib.loads((ROOT / "uv.lock").read_text(encoding="utf-8"))
    (locked_project,) = (
        package
        for package in lock["package"]
        if package["name"] == project["name"] and package["source"] == {"editable": "."}
    )
    assert locked_project["version"] == expected


def test_citation_has_top_level_authors_and_standard_license() -> None:
    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    preferred_citation_offset = citation.index("preferred-citation:")
    assert "\nauthors:\n" in citation[:preferred_citation_offset]
    assert (ROOT / "LICENSE").is_file()
    assert not (ROOT / "MIT License").exists()

import re
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
    assert _declared_version(ROOT / "pyproject.toml") == "0.1.4"
    assert _declared_version(ROOT / "CITATION.cff") == "0.1.4"
    assert __version__ == "0.1.4"


def test_citation_has_top_level_authors_and_standard_license() -> None:
    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    preferred_citation_offset = citation.index("preferred-citation:")
    assert "\nauthors:\n" in citation[:preferred_citation_offset]
    assert (ROOT / "LICENSE").is_file()
    assert not (ROOT / "MIT License").exists()

import csv
import json
import re
import shutil
from pathlib import Path
from urllib.parse import unquote, urlsplit

import pytest
from documentation import (
    ANNOTATION_TUTORIAL,
    DOCUMENTS,
    MARKDOWN,
    MINIMAL_TUTORIAL,
    ROOT,
    cli_commands,
    example_arguments,
    generate_annotation_inputs,
)

from kffractbias.cli import build_parser, main


@pytest.mark.parametrize("document", DOCUMENTS, ids=lambda path: str(path.relative_to(ROOT)))
def test_documented_cli_commands_parse(document, monkeypatch):
    # Most README paths are placeholders; syntax checks must not require real genomes.
    monkeypatch.setattr("kffractbias.cli._path", Path)
    parser = build_parser()
    for command in cli_commands(document):
        try:
            parser.parse_args(command.argv)
        except SystemExit as exc:
            pytest.fail(f"Invalid CLI example at {document}:{command.line}: {exc}")


def heading_ids(path):
    identifiers = set()
    tokens = MARKDOWN.parse(path.read_text(encoding="utf-8"))
    for index, token in enumerate(tokens):
        if token.type != "heading_open":
            continue
        text = "".join(
            child.content
            for child in tokens[index + 1].children or []
            if child.type in {"text", "code_inline"}
        )
        base = re.sub(r"[^\w\s-]", "", text.lower()).replace(" ", "-")
        identifier, suffix = base, 0
        while identifier in identifiers:
            suffix += 1
            identifier = f"{base}-{suffix}"
        identifiers.add(identifier)
    return identifiers


@pytest.mark.parametrize("document", DOCUMENTS, ids=lambda path: str(path.relative_to(ROOT)))
def test_internal_document_links_exist(document):
    for token in MARKDOWN.parse(document.read_text(encoding="utf-8")):
        for child in token.children or []:
            if child.type not in {"link_open", "image"}:
                continue
            href = child.attrGet("href" if child.type == "link_open" else "src")
            url = urlsplit(href)
            if url.scheme or url.netloc:
                continue
            target = (document.parent / unquote(url.path)).resolve() if url.path else document
            assert target.is_relative_to(ROOT), f"Link escapes repository: {document}: {href}"
            assert target.exists(), f"Broken link: {document}: {href}"
            if url.fragment and target.suffix == ".md":
                assert unquote(url.fragment) in heading_ids(target), (
                    f"Broken heading link: {document}: {href}"
                )


def test_minimal_documented_command_produces_expected_retention(tmp_path, monkeypatch):
    shutil.copytree(ROOT / "examples" / "minimal", tmp_path / "examples" / "minimal")
    monkeypatch.chdir(tmp_path)
    arguments = example_arguments(MINIMAL_TUTORIAL, "calculate")
    options = build_parser().parse_args(arguments)
    assert main(arguments) == 0
    output = options.output_dir
    summary = json.loads((output / f"{options.prefix}.summary.json").read_text())
    assert summary["parameters"]["denominator"] == "all"
    assert summary["parameters"]["synteny_format"] == "jcvi"
    assert summary["counts"]["synteny_pair_count"] == 3
    assert summary["counts"]["gene_table_row_count"] == 8
    assert summary["counts"]["window_table_row_count"] == 6
    with (output / f"{options.prefix}.windows.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {
        (row["query_seqid"], int(row["start_rank"])): float(row["retention_percent"])
        for row in rows
    } == {
        ("chrA", 1): 100,
        ("chrA", 2): 50,
        ("chrA", 3): 0,
        ("chrB", 1): 0,
        ("chrB", 2): 50,
        ("chrB", 3): 50,
    }
    assert not list(output.glob("*.plot.*"))


def test_annotation_generator_and_documented_validation(tmp_path, monkeypatch, capsys):
    generate_annotation_inputs(tmp_path)
    monkeypatch.chdir(tmp_path)
    assert main(example_arguments(ANNOTATION_TUTORIAL, "validate")) == 0
    result = json.loads(capsys.readouterr().out)
    for label, count in (("target", 8), ("query", 16)):
        assert result[label]["matched_gene_count"] == count
        assert result[label]["selected_gene_count"] == count
        assert result[label]["collapsed_isoform_count"] == 0
        assert result[label]["unresolved_locus_count"] == 0

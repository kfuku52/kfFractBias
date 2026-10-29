import csv
import json
from pathlib import Path

import pytest

from kffractbias.cli import main
from kffractbias.comparison import (
    ComparisonInput,
    aligned_windows,
    build_loci,
    comparison_report,
)


def fixture(tmp_path: Path, name="first", missing=False):
    target = tmp_path / "reference.bed"
    target.write_text(
        "".join(f"chr{c}\t{i * 100}\t{i * 100 + 90}\tt{c}_{i}\n" for c in [1, 2] for i in range(6))
    )
    query = tmp_path / (name + ".bed")
    mapping = tmp_path / (name + ".mapping.tsv")
    anchors = tmp_path / (name + ".anchors")
    beds, maps, blocks = [], [], []
    for c in [1, 2]:
        for copy, indices in [("A", range(6)), ("B", [0, 1, 3, 5])]:
            fragment = f"{name}_{c}_{copy}"
            pairs = []
            for i in indices:
                if missing and c == 2 and copy == "B":
                    continue
                gene = f"{name}{c}{copy}{i}"
                beds.append(f"{fragment}\t{i * 100}\t{i * 100 + 90}\t{gene}\n")
                maps.append(f"{gene}\tG{c}\t{copy}\t{fragment}\n")
                pairs.append(f"t{c}_{i}\t{gene}\t100\n")
            blocks.append("###\n" + "".join(pairs))
    query.write_text("".join(beds))
    mapping.write_text("gene_id\tgroup_id\tcopy\tsegment_id\n" + "".join(maps))
    anchors.write_text("".join(blocks))
    return ComparisonInput(name, name, 2, target, query, anchors, mapping)


def manifest(tmp_path, items):
    path = tmp_path / "manifest.tsv"
    with path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "analysis_id",
                "species",
                "expected_copies",
                "target_bed",
                "query_bed",
                "synteny_file",
                "mapping_file",
            ],
            delimiter="\t",
        )
        w.writeheader()
        for item in items:
            w.writerow(vars(item))
    return path


def test_counts_zero_retention_inside_a_callable_interior(tmp_path):
    item = fixture(tmp_path)
    _, loci, audit = build_loci(item, minimum_anchors=2, trim=0)
    assert [r["retained_copies"] for r in loci[:6]] == [2, 2, 1, 2, 1, 2]
    assert all(r["callable"] for r in loci)
    assert audit["callable_interiors"] == 4
    windows, shared = aligned_windows(
        [loci], window_size=3, step_size=3, denominator="all", minimum_callable_fraction=1
    )
    assert len(windows[0]) == 4  # Never bridge reference chromosomes.
    assert windows[0][0]["total_percent"] == pytest.approx(500 / 3)
    assert windows[0][0]["normalised_percent"] == pytest.approx(250 / 3)
    assert windows[0][0]["copy_percent"] == {"G1:A": 100, "G1:B": pytest.approx(200 / 3)}
    assert all(shared)


def test_duplicate_fragment_hits_are_one_copy(tmp_path):
    item = fixture(tmp_path)
    item.query_bed.write_text(item.query_bed.read_text() + "extra_fragment\t0\t90\textra\n")
    item.mapping_file.write_text(item.mapping_file.read_text() + "extra\tG1\tA\textra_fragment\n")
    item.synteny_file.write_text(
        item.synteny_file.read_text() + "###\nt1_0\textra\t100\nt1_0\textra\t100\n"
    )
    _, loci, _ = build_loci(item, minimum_anchors=2, trim=0)
    assert loci[0]["retained_copies"] == 2
    assert loci[0]["copies"]["G1:A"] == 1


def test_missing_copy_is_unknown_and_shared_windows_keep_one_denominator(tmp_path):
    a = fixture(tmp_path, "a")
    b = fixture(tmp_path, "b")
    # Keep the independent copy assignment but remove its collinear evidence.
    b.synteny_file.write_text(
        b.synteny_file.read_text().replace("t2_", "t1_").replace("\tb2", "\tb2")
    )
    _, pa, _ = build_loci(a, minimum_anchors=2, trim=0)
    _, pb, _ = build_loci(b, minimum_anchors=2, trim=0)
    windows, shared = aligned_windows(
        [pa, pb], window_size=3, step_size=3, denominator="all", minimum_callable_fraction=1
    )
    assert shared[6:] == [False] * 6
    assert all(r["total_percent"] is None for p in windows for r in p[2:])
    assert [r["common_callable_genes"] for r in windows[0]] == [
        r["common_callable_genes"] for r in windows[1]
    ]


def test_ambiguous_local_groups_are_not_pooled(tmp_path):
    item = fixture(tmp_path)
    item.synteny_file.write_text(
        item.synteny_file.read_text()
        + "###\n"
        + "".join(
            f"t1_{i}\tfirst2{copy}{i}\t100\n"
            for copy, indices in [("A", range(6)), ("B", [0, 1, 3, 5])]
            for i in indices
        )
    )
    _, loci, audit = build_loci(item, minimum_anchors=2, trim=0)
    assert all(not r["callable"] for r in loci[:6])
    assert audit["ambiguous_loci"] == 6


def test_syntenic_denominator_is_intersection_and_does_not_cross_chromosomes(tmp_path):
    item = fixture(tmp_path)
    _, p, _ = build_loci(item, minimum_anchors=2, trim=0)
    q = [
        {**r, "retained_copies": 0 if r["target_gene"] == "t1_1" else r["retained_copies"]}
        for r in p
    ]
    windows, _ = aligned_windows(
        [p, q], window_size=3, step_size=1, denominator="syntenic", minimum_callable_fraction=1
    )
    assert len(windows[0]) == 7
    assert windows[0][0]["start_gene"] == "t1_0"
    assert windows[0][0]["end_gene"] == "t1_3"


def test_report_cli_and_hashed_tables(tmp_path):
    item = fixture(tmp_path)
    path = manifest(tmp_path, [item])
    out = tmp_path / "output"
    assert (
        main(
            [
                "report",
                "--manifest",
                str(path),
                "--output-dir",
                str(out),
                "--window-size",
                "3",
                "--minimum-anchors",
                "2",
                "--trim",
                "0",
                "--no-plot",
            ]
        )
        == 0
    )
    summary = json.loads((out / "comparison.report.json").read_text())
    assert summary["estimable_windows_per_analysis"] == 8
    assert summary["common_callable_loci"] == 12
    assert set(summary["output_sha256"]) == {"loci", "windows"}
    assert len(summary["inputs"]) == 5


@pytest.mark.parametrize(
    "problem", ["unknown_gene", "missing_copy", "reference_coordinates", "ragged", "collision"]
)
def test_invalid_inputs_preserve_existing_report(tmp_path, problem):
    a = fixture(tmp_path, "a")
    b = fixture(tmp_path, "b")
    path = manifest(tmp_path, [a, b])
    out = tmp_path / "output"
    comparison_report(path, out, window_size=3, minimum_anchors=2, trim=0)
    previous = (out / "comparison.report.json").read_bytes()
    if problem == "unknown_gene":
        a.synteny_file.write_text("unknown\ta1A0\t10\n")
    elif problem == "missing_copy":
        a.mapping_file.write_text(a.mapping_file.read_text().replace("\tB\t", "\tA\t"))
    elif problem == "reference_coordinates":
        target = tmp_path / "different.bed"
        target.write_text(a.target_bed.read_text().replace("\t90\t", "\t91\t"))
        path.write_text(path.read_text().replace(str(b.target_bed), str(target), 1))
    elif problem == "ragged":
        path.write_text(path.read_text() + "extra\trow\n")
    else:
        path = out / "comparison.loci.tsv"
    with pytest.raises((ValueError, OSError)):
        comparison_report(path, out, window_size=3, minimum_anchors=2, trim=0)
    assert (out / "comparison.report.json").read_bytes() == previous


def test_plots_and_all_reference_pages(tmp_path):
    pytest.importorskip("matplotlib")
    a = fixture(tmp_path)
    path = manifest(tmp_path, [a])
    summary = comparison_report(
        path,
        tmp_path / "output",
        window_size=3,
        minimum_anchors=2,
        trim=0,
        plot_config={
            "representative_seqids": ["chr1", "chr2"],
            "reference_seqids": ["chr2", "chr1"],
            "width_pt": 650,
            "height_pt": 260,
        },
    )
    assert summary["plot"]["all_reference_pages"] == 2
    assert summary["plot"]["font_size_pt"] == 8
    assert len(summary["output_sha256"]) == 9


def test_callable_windows_select_shared_opportunities_before_ranking(tmp_path):
    item = fixture(tmp_path)
    _, p, _ = build_loci(item, minimum_anchors=2, trim=0)
    q = [{**r, "callable": 0 if r["target_gene"] == "t1_1" else r["callable"]} for r in p]
    windows, _ = aligned_windows(
        [p, q], window_size=3, step_size=1, denominator="callable", minimum_callable_fraction=1
    )
    assert len(windows[0]) == 7
    assert windows[0][0]["common_callable_genes"] == 3
    assert windows[0][0]["start_rank"] == 1 and windows[0][0]["end_rank"] == 4
    assert [r["start_gene"] for r in windows[0]] == [r["start_gene"] for r in windows[1]]


def test_manifest_change_between_discovery_and_snapshot_is_rejected(tmp_path, monkeypatch):
    import kffractbias.comparison as module

    item = fixture(tmp_path)
    path = manifest(tmp_path, [item])
    original = module.read_manifest

    def change_after_discovery(source):
        items = original(source)
        source.write_text(source.read_text().replace("first\tfirst", "other\tother"))
        return items

    monkeypatch.setattr(module, "read_manifest", change_after_discovery)
    with pytest.raises(RuntimeError, match="Manifest changed"):
        comparison_report(path, tmp_path / "out", window_size=3, minimum_anchors=2, trim=0)
    assert not (tmp_path / "out/comparison.report.json").exists()


def test_plot_configuration_file_is_captured_and_validated(tmp_path):
    item = fixture(tmp_path)
    path = manifest(tmp_path, [item])
    config = tmp_path / "plot.json"
    config.write_text("[]")
    with pytest.raises(ValueError, match="object"):
        comparison_report(path, tmp_path / "out", plot_config_file=config)
    config.write_text('{"representative_seqids":["chr1"],"width_pt":650,"height_pt":260}')
    report = comparison_report(
        path, tmp_path / "out", window_size=3, minimum_anchors=2, trim=0, plot_config_file=config
    )
    assert report["inputs"]["plot_config"]["path"] == str(config)
    assert report["plot"]["representative_seqids"] == ["chr1"]


def test_five_copies_are_normalised_by_five_with_zero_and_duplicate_hits(tmp_path):
    from dataclasses import replace

    item = fixture(tmp_path)
    original_bed = item.query_bed.read_text()
    original_map = item.mapping_file.read_text().splitlines()[1:]
    original_blocks = item.synteny_file.read_text().split("###\n")[1:]
    for copy in "CDE":
        item.query_bed.write_text(
            item.query_bed.read_text()
            + "".join(
                line.replace("_A", "_" + copy).replace("A", copy)
                for line in original_bed.splitlines(keepends=True)
                if "_A\t" in line
            )
        )
        item.mapping_file.write_text(
            item.mapping_file.read_text()
            + "".join(
                line.replace("A", copy).replace("_A", "_" + copy) + "\n"
                for line in original_map
                if "\tA\t" in line
            )
        )
        item.synteny_file.write_text(
            item.synteny_file.read_text()
            + "".join(
                "###\n" + block.replace("A", copy)
                for block in original_blocks
                if "\tfirst1A" in block or "\tfirst2A" in block
            )
        )
    item.synteny_file.write_text(item.synteny_file.read_text() + "###\nt1_0\tfirst1E0\t100\n")
    _, p, _ = build_loci(replace(item, expected_copies=5), minimum_anchors=2, trim=0)
    assert [r["retained_copies"] for r in p[:6]] == [5, 5, 4, 5, 4, 5]
    windows, _ = aligned_windows(
        [p], window_size=3, step_size=3, denominator="all", minimum_callable_fraction=1
    )
    assert windows[0][0]["normalised_percent"] == pytest.approx(100 * 14 / 15)
    assert windows[0][0]["total_percent"] == pytest.approx(100 * 14 / 3)

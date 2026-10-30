"""Publication styling must preserve sequence identity and numerical profiles."""

import json

import pytest

from kffractbias.comparison import ComparisonInput
from kffractbias.comparison_plot import plot_comparison


def inputs(tmp_path, series_count=12):
    reference = tmp_path / "reference.bed"
    query = tmp_path / "query.bed"
    reference.write_text(
        "".join(f"chr{c}\t{i * 100}\t{i * 100 + 90}\tt{c}_{i}\n" for c in [1, 2] for i in range(6))
    )
    query.write_text("".join(f"q{i}\t0\t10\tg{i}\n" for i in range(1, series_count + 2)))
    item = ComparisonInput("a", "Example_species", 2, reference, query, tmp_path / "anchors", None)
    rows = []
    for chrom, first in [("chr1", 1), ("chr2", 2)]:
        for start in [1, 4]:
            values = {
                f"all:q{i}": (60 if i == 2 else 40 if i == 3 else 1)
                for i in range(first, series_count + 1)
            }
            rows.append(
                {
                    "target_seqid": chrom,
                    "start_rank": start,
                    "end_rank": start + 2,
                    "total_percent": sum(values.values()),
                    "normalised_percent": sum(values.values()) / 2,
                    "copy_percent": values,
                    "segment_ids": {k: k for k in values},
                }
            )
    return item, rows


def draw(tmp_path, config, series_count=12):
    item, rows = inputs(tmp_path, series_count)
    out = tmp_path / "plots"
    out.mkdir()
    outputs = {
        f"{kind}_{ext}": out / f"{kind}.{ext}"
        for kind in ["normalised", "raw"]
        for ext in ["pdf", "png", "svg"]
    }
    outputs["all_pdf"] = out / "all.pdf"
    return plot_comparison(
        [item],
        [rows],
        {p: p for p in outputs.values()},
        outputs,
        {"width_pt": 650, "height_pt": 300, **config},
    )


@pytest.mark.parametrize("palette", ["categorical", "pastel"])
def test_colours_are_unique_above_ten_and_stable_when_a_panel_lacks_a_series(
    tmp_path, monkeypatch, palette
):
    import matplotlib.figure
    from matplotlib.colors import to_hex

    snapshots = []
    original = matplotlib.figure.Figure.savefig

    def capture(figure, *args, **kwargs):
        snapshots.extend(
            (
                axis.get_title(),
                {
                    line.get_label(): to_hex(line.get_color())
                    for line in axis.lines
                    if line.get_label().startswith("all:")
                },
            )
            for axis in figure.axes
            if axis.get_xlabel() == "Gene rank"
        )
        return original(figure, *args, **kwargs)

    monkeypatch.setattr(matplotlib.figure.Figure, "savefig", capture)
    result = draw(tmp_path, {"palette": palette}, series_count=40)
    colours = result["series_colours"]["a"]
    assert len(colours) == len(set(colours.values())) == 40
    assert "q41" not in colours  # No evidence, but its input label remains valid for an override.
    for _, seen in snapshots:
        for key, value in seen.items():
            assert value == colours[key.split(":")[1]]
    assert any(title == "chr1" and "all:q1" in seen for title, seen in snapshots)
    assert any(title == "chr2" and "all:q1" not in seen for title, seen in snapshots)


@pytest.mark.parametrize("scale,total", [("raw", 110), ("normalised", 55)])
def test_graphical_legend_overrides_and_all_reference_scale(tmp_path, monkeypatch, scale, total):
    import matplotlib.backends.backend_pdf
    from matplotlib.colors import to_hex

    pages = []
    original = matplotlib.backends.backend_pdf.PdfPages.savefig

    def capture(pdf, figure, *args, **kwargs):
        pages.append(
            (
                figure.axes[0].lines[0].get_ydata().tolist(),
                [t.get_text() for t in figure.legends[0].get_texts()],
                [
                    (to_hex(line.get_color()), line.get_alpha())
                    for line in figure.axes[0].lines
                    if line.get_label() == "all:q2"
                ],
            )
        )
        return original(pdf, figure, *args, **kwargs)

    monkeypatch.setattr(matplotlib.backends.backend_pdf.PdfPages, "savefig", capture)
    result = draw(
        tmp_path,
        {
            "palette": "pastel",
            "series_colours": {"a": {"q2": "#8bb9dc", "q3": "#e6b979", "q13": "#e6b979"}},
            "total_colour": "#b53c31",
            "series_alpha": 0.6,
            "legend": "graphical",
            "all_profiles_scale": scale,
        },
    )
    assert pages[0][0] == [total, total]
    assert pages[0][1] == [
        "Total" if scale == "raw" else "Total / expected copies",
        "Chromosome/contig",
        "100% reference",
    ]
    assert all(page[2] == [("#8bb9dc", 0.6)] for page in pages)
    assert result["series_colours"]["a"]["q2"] == to_hex("#8bb9dc")
    assert result["unused_series_colour_labels"] == {"a": ["q13"]}
    assert result["total_colour"] == "#b53c31"
    assert result["series_alpha"] == 0.6
    assert result["all_profiles_scale"] == scale


@pytest.mark.parametrize(
    "bad",
    [
        {"palette": "unknown"},
        {"legend": "unknown"},
        {"all_profiles_scale": "unknown"},
        {"series_alpha": 0},
        {"series_alpha": float("nan")},
        {"series_alpha": True},
        {"total_colour": "bad colour"},
        {"total_colour": "none"},
        {"series_colours": []},
        {"series_colours": {"missing_analysis": {}}},
        {"series_colours": {"a": {"missing_sequence": "red"}}},
        {"series_colours": {"a": {"q1": "#ff000080"}}},
    ],
)
def test_invalid_style_is_rejected(tmp_path, bad):
    with pytest.raises(ValueError):
        draw(tmp_path, bad)
    assert not list((tmp_path / "plots").iterdir())


def test_palette_is_anchored_to_input_labels_when_evidence_is_filtered(tmp_path):
    item, rows = inputs(tmp_path)
    config = {"width_pt": 650, "height_pt": 300, "legend": "none"}
    outputs = {
        f"{k}_{e}": tmp_path / f"{k}.{e}"
        for k in ["normalised", "raw"]
        for e in ["pdf", "png", "svg"]
    }
    outputs["all_pdf"] = tmp_path / "all.pdf"
    paths = {p: p for p in outputs.values()}
    a = plot_comparison([item], [rows], paths, outputs, config)
    b = plot_comparison(
        [item], [[r for r in rows if r["target_seqid"] == "chr2"]], paths, outputs, config
    )
    for label, value in b["series_colours"]["a"].items():
        assert value == a["series_colours"]["a"][label]
    assert b["all_profiles_scale"] == "normalised"  # Backwards-compatible default.
    assert json.loads(json.dumps(b))["legend"] == "none"

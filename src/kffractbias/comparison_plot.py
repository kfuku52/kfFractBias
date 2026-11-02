"""Publication plots for reference-aligned copy profiles; no inference here."""

from __future__ import annotations

import colorsys
import math
from collections import Counter
from pathlib import Path
from typing import Any

from .comparison import ComparisonInput, table
from .io import natural_key, read_bed, sha256_file


def plot_comparison(
    items: list[ComparisonInput],
    windows: list[list[dict[str, Any]]],
    staged: dict[Path, Path | None],
    outputs: dict[str, Path],
    config: dict[str, Any],
) -> dict[str, Any]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import to_hex, to_rgba
    from matplotlib.legend_handler import HandlerTuple
    from matplotlib.lines import Line2D
    from matplotlib.text import Text
    from matplotlib.ticker import FuncFormatter

    allowed = {
        "font_family",
        "font_size",
        "font_files",
        "width_pt",
        "height_pt",
        "dpi",
        "representative_seqids",
        "reference_seqids",
        "palette",
        "series_colours",
        "total_colour",
        "series_alpha",
        "legend",
        "all_profiles_scale",
    }
    if set(config) - allowed:
        raise ValueError("Unsupported comparison plot configuration")
    palette = config.get("palette", "categorical")
    legend = config.get("legend", "text")
    all_scale = config.get("all_profiles_scale", "normalised")
    if palette not in ("categorical", "pastel"):
        raise ValueError("palette must be categorical or pastel")
    if legend not in ("text", "graphical", "none"):
        raise ValueError("legend must be text, graphical or none")
    if all_scale not in ("normalised", "raw"):
        raise ValueError("all_profiles_scale must be normalised or raw")
    alpha = config.get("series_alpha", 0.85)
    if isinstance(alpha, bool) or not isinstance(alpha, (int, float)) or not 0 < alpha <= 1:
        raise ValueError("series_alpha must be a finite number in (0, 1]")

    def colour(value: Any) -> str:
        if not isinstance(value, str):
            raise ValueError("Colours must be opaque Matplotlib colour strings")
        try:
            rgba = to_rgba(value)
        except ValueError as exc:
            raise ValueError(f"Invalid colour: {value}") from exc
        if rgba[3] != 1:
            raise ValueError("Colours must be opaque; use series_alpha for transparency")
        return str(to_hex(rgba))

    total_colour = colour(config.get("total_colour", "black"))
    overrides = config.get("series_colours", {})
    if not isinstance(overrides, dict) or set(overrides) - {i.analysis_id for i in items}:
        raise ValueError("series_colours must map known analysis IDs to colour objects")

    def default_colour(index: int) -> str:
        # Preserve the first ten categorical colours, then extend without cycling.
        rgb = (
            plt.get_cmap("tab10")(index)[:3]
            if index < 10
            else colorsys.hls_to_rgb((index * 0.618033988749895) % 1, 0.5, 0.7)
        )
        if palette == "pastel":
            rgb = tuple(0.55 + 0.45 * c for c in rgb)
        return str(to_hex(rgb))

    series_colours: dict[str, dict[str, str]] = {}
    unused_colours: dict[str, list[str]] = {}
    for item, rows in zip(items, windows, strict=True):
        # Use the complete input label universe, not the labels on one reference
        # chromosome. Filtering windows or adding a panel cannot shift a colour.
        labels = (
            {g.seqid for g in read_bed(item.query_bed)}
            if item.mapping_file is None
            else {
                r["copy"]
                for r in table(item.mapping_file, {"gene_id", "group_id", "copy", "segment_id"})
            }
        )
        requested = overrides.get(item.analysis_id, {})
        if not isinstance(requested, dict) or set(requested) - labels:
            raise ValueError("series_colours contains an unknown copy/sequence label")
        observed_labels = {k.split(":", 1)[1] for r in rows for k in r["copy_percent"]}
        if observed_labels - labels:
            raise ValueError("Plotted copy labels are absent from comparison inputs")
        series_colours[item.analysis_id] = {
            label: colour(requested[label]) if label in requested else default_colour(index)
            for index, label in enumerate(sorted(labels, key=natural_key))
            if label in observed_labels or label in requested
        }
        unused_colours[item.analysis_id] = sorted(set(requested) - observed_labels, key=natural_key)
    family = config.get("font_family", "DejaVu Sans")
    size = float(config.get("font_size", 8))
    width, height = (
        float(config.get("width_pt", 518.4)),
        float(config.get("height_pt", 180 * len(items))),
    )
    dpi = int(config.get("dpi", 300))
    if not all(math.isfinite(v) and v > 0 for v in [size, width, height, dpi]):
        raise ValueError("Plot dimensions, DPI and font size must be positive and finite")
    selector = font_manager.FontManager()
    for path in config.get("font_files", []):
        selector.addfont(str(path))
    selector.findfont(font_manager.FontProperties(family=family), fallback_to_default=False)
    observed = {r["target_seqid"] for profile in windows for r in profile}
    reference_counts = Counter(g.seqid for g in read_bed(items[0].target_bed))
    reference_names = set(reference_counts)
    chromosomes = config.get(
        "reference_seqids", sorted(observed or reference_names, key=natural_key)
    )
    if (
        not chromosomes
        or len(set(chromosomes)) != len(chromosomes)
        or observed - set(chromosomes)
        or set(chromosomes) - reference_names
    ):
        raise ValueError("Unknown, duplicate or incomplete reference sequence display order")
    representatives = config.get("representative_seqids", chromosomes[:3])
    if (
        not representatives
        or len(set(representatives)) != len(representatives)
        or set(representatives) - set(chromosomes)
    ):
        raise ValueError("Unknown or duplicate representative reference sequences")
    fonts: dict[str, str] = {}

    def finish(figure: Any) -> None:
        for text in figure.findobj(Text):
            prop = text.get_fontproperties().copy()
            prop.set_size(size)
            prop.set_family(family)
            prop.set_file(selector.findfont(prop, fallback_to_default=False))
            text.set_fontproperties(prop)
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        unused = set()
        for axis in figure.axes:
            for ticks, limits in [(axis.xaxis, axis.get_xlim()), (axis.yaxis, axis.get_ylim())]:
                low, high = sorted(limits)
                for tick in ticks.get_major_ticks() + ticks.get_minor_ticks():
                    if not low <= tick.get_loc() <= high:
                        unused.update([tick.label1, tick.label2])
        for text in figure.findobj(Text):
            if not text.get_visible() or not text.get_text() or text in unused:
                continue
            path = text.get_fontproperties().get_file()
            if path is None or font_manager.FontProperties(fname=path).get_name() != family:
                raise ValueError("Figure resolved an unexpected font")
            fonts[path] = sha256_file(path)
            box = text.get_window_extent(renderer)
            if (
                box.x0 < -1
                or box.y0 < -1
                or box.x1 > figure.bbox.width + 1
                or box.y1 > figure.bbox.height + 1
            ):
                raise ValueError(
                    f"Text outside figure: {text.get_text()}; increase plot dimensions"
                )

    def curves(
        axis: Any, item: ComparisonInput, rows: list[dict[str, Any]], chrom: str, normalised: bool
    ) -> None:
        selected = [r for r in rows if r["target_seqid"] == chrom]
        x = [(r["start_rank"] + r["end_rank"]) / 2 for r in selected]
        field = "normalised_percent" if normalised else "total_percent"
        y = [r[field] if r[field] is not None else math.nan for r in selected]
        # The numerical rank remains fixed; NA windows leave visible gaps.
        axis.plot(
            x,
            y,
            color=total_colour,
            linewidth=1.1,
            zorder=3,
            label="Total / expected copies" if normalised else "Total",
        )
        keys = sorted({k for r in selected for k in r["copy_percent"]}, key=natural_key)
        for key in keys:
            segment_values = {r["segment_ids"][key] for r in selected if key in r["segment_ids"]}
            for segment in sorted(segment_values):
                values = [
                    r["copy_percent"].get(key, math.nan)
                    if r["segment_ids"].get(key) == segment
                    else math.nan
                    for r in selected
                ]
                axis.plot(
                    x,
                    values,
                    color=series_colours[item.analysis_id][key.split(":", 1)[1]],
                    linewidth=0.55,
                    alpha=alpha,
                    label=key,
                )
        axis.axhline(100, color=".6", linestyle=":", linewidth=0.5)
        axis.set_title(chrom)
        axis.set_xlabel("Gene rank")
        count = reference_counts[chrom]
        axis.set_xlim((1, count) if count > 1 else (0.5, 1.5))
        axis.xaxis.set_major_locator(plt.MaxNLocator(4, integer=True))
        axis.xaxis.set_major_formatter(
            FuncFormatter(
                lambda value, position: f"{value / 1000:g}k" if abs(value) >= 1000 else f"{value:g}"
            )
        )
        axis.spines[["top", "right"]].set_visible(False)
        if not any(math.isfinite(v) for v in y):
            axis.text(0.5, 0.5, "Not estimable", transform=axis.transAxes, ha="center")
        upper = normalised_upper if normalised else raw_upper
        axis.set_ylim(0, upper)

    def add_legend(figure: Any, normalised: bool) -> None:
        if legend == "text":
            name = "Black" if total_colour == "#000000" else "Total"
            figure.supxlabel(
                f"{name}: sum / expected depth; colours: local copies (unphased)"
                if normalised
                else f"{name}: summed copies; colours: local copies (unphased)"
            )
        elif legend == "graphical":
            samples = list(dict.fromkeys(c for m in series_colours.values() for c in m.values()))[
                :3
            ]
            if not samples:
                samples = [default_colour(0)]
            handles = [
                Line2D([], [], color=total_colour, linewidth=1.1),
                tuple(Line2D([], [], color=c, alpha=alpha, linewidth=0.8) for c in samples),
                Line2D([], [], color=".6", linestyle=":", linewidth=0.5),
            ]
            figure.legend(
                handles=handles,
                labels=[
                    "Total / expected copies" if normalised else "Total",
                    "Chromosome/contig"
                    if all(i.mapping_file is None for i in items)
                    else "Local copy",
                    "100% reference",
                ],
                handler_map={tuple: HandlerTuple(ndivide=None)},
                loc="outside lower center",
                ncol=3,
                frameon=False,
            )

    style = {
        "font.family": family,
        "font.size": size,
        "axes.labelsize": size,
        "axes.titlesize": size,
        "xtick.labelsize": size,
        "ytick.labelsize": size,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
    with plt.rc_context(style):
        normalised_upper = max(
            102,
            max(
                (
                    r["normalised_percent"]
                    for p in windows
                    for r in p
                    if r["normalised_percent"] is not None
                ),
                default=100,
            )
            * 1.1,
        )
        raw_upper = max(
            102,
            max(
                (r["total_percent"] for p in windows for r in p if r["total_percent"] is not None),
                default=100,
            )
            * 1.1,
        )
        for normalised in [True, False]:
            kind = "normalised" if normalised else "raw"
            figure = plt.figure(figsize=(width / 72, height / 72), layout="constrained")
            try:
                grid = figure.add_gridspec(len(items), 2, width_ratios=[1.15, 1])
                for index, (item, rows) in enumerate(zip(items, windows, strict=True)):
                    axis = figure.add_subplot(grid[index, 0])
                    field = "normalised_percent" if normalised else "total_percent"
                    values = [
                        [
                            r[field]
                            for r in rows
                            if r["target_seqid"] == chrom and r[field] is not None
                        ]
                        for chrom in chromosomes
                    ]
                    axis.boxplot(
                        [v or [math.nan] for v in values],
                        tick_labels=chromosomes,
                        showfliers=False,
                        widths=0.6,
                        medianprops={"color": "black"},
                        patch_artist=True,
                        boxprops={"facecolor": "white"},
                    )
                    axis.tick_params(axis="x", labelrotation=90)
                    axis.set_title(item.species.replace("_", " "), fontstyle="italic", loc="left")
                    axis.set_ylabel(
                        "Expected-copy detection (%)"
                        if normalised
                        else "Summed syntelog detection (%)"
                    )
                    axis.set_xlabel("Reference chromosome")
                    axis.spines[["top", "right"]].set_visible(False)
                    if normalised:
                        axis.set_ylim(0, normalised_upper)
                    else:
                        axis.set_ylim(0, raw_upper)
                    if not any(values):
                        axis.text(0.5, 0.5, "Not estimable", transform=axis.transAxes, ha="center")
                    sub = grid[index, 1].subgridspec(1, len(representatives))
                    for j, chrom in enumerate(representatives):
                        trace = figure.add_subplot(sub[0, j])
                        curves(trace, item, rows, chrom, normalised)
                        if j == 0:
                            trace.set_ylabel("Syntelog detection (%)")
                        trace.xaxis.set_major_locator(plt.MaxNLocator(2))
                add_legend(figure, normalised)
                finish(figure)
                for ext in ["pdf", "png", "svg"]:
                    path = staged[outputs[kind + "_" + ext]]
                    if path is None:
                        raise ValueError("Missing staged figure path")
                    figure.savefig(path, dpi=dpi)
            finally:
                plt.close(figure)
        all_path = staged[outputs["all_pdf"]]
        if all_path is None:
            raise ValueError("Missing all-chromosome PDF path")
        with PdfPages(all_path, metadata={"CreationDate": None, "ModDate": None}) as pdf:
            for chrom in chromosomes:
                figure, axes = plt.subplots(
                    len(items),
                    1,
                    figsize=(width / 72, height / 72),
                    squeeze=False,
                    layout="constrained",
                )
                try:
                    for item, rows, axis in zip(items, windows, axes[:, 0], strict=True):
                        curves(axis, item, rows, chrom, all_scale == "normalised")
                        axis.set_title(
                            item.species.replace("_", " ") + ": " + chrom,
                            fontstyle="italic",
                            loc="left",
                        )
                        axis.set_ylabel("Syntelog detection (%)")
                    add_legend(figure, all_scale == "normalised")
                    finish(figure)
                    pdf.savefig(figure)
                finally:
                    plt.close(figure)
    return {
        "font_family": family,
        "font_size_pt": size,
        "font_hashes": fonts,
        "figure_size_pt": [width, height],
        "representative_seqids": representatives,
        "reference_seqids": chromosomes,
        "all_reference_pages": len(chromosomes),
        "copies_are_local": True,
        "reference_gene_limits": {chrom: [1, reference_counts[chrom]] for chrom in chromosomes},
        "rank_tick_labels": "k denotes 1000 reference genes; axes retain full reference ranks",
        "total_black": "normalised total divides the copy count by expected_copies",
        "coloured_tracks": "individual local copies; lines break at native segment changes and missing windows",
        "palette": palette,
        "series_colours": series_colours,
        "unused_series_colour_labels": unused_colours,
        "total_colour": total_colour,
        "series_alpha": alpha,
        "legend": legend,
        "all_profiles_scale": all_scale,
    }

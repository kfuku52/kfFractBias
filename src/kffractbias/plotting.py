from __future__ import annotations

import math
from array import array
from collections.abc import Iterable, Iterator
from itertools import groupby, islice
from pathlib import Path

from .profiles import WindowRow

MAX_PANELS_PER_PAGE = 6
MAX_QUERIES_PER_PANEL = 12


def preflight_plot() -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "Plotting requires matplotlib; rerun with --no-plot to omit figures"
        ) from exc


def _panels(rows: Iterable[WindowRow]) -> Iterator[tuple[str, list[tuple[str, array, array]]]]:
    # RetentionProfile emits target, query, window order. Retain numeric arrays
    # for one bounded page, never a copy of the entire window table.
    for target, target_rows in groupby(rows, key=lambda row: row["target_seqid"]):
        series: list[tuple[str, array, array]] = []
        for query, query_rows in groupby(target_rows, key=lambda row: row["query_seqid"]):
            x, y = array("d"), array("d")
            for row in query_rows:
                x.append((row["start_rank"] + row["end_rank"]) / 2)
                y.append(float(row["retention_percent"]))
            series.append((query, x, y))
            if len(series) == MAX_QUERIES_PER_PANEL:
                yield target, series
                series = []
        if series:
            yield target, series


def plot_windows(
    rows: Iterable[WindowRow],
    pdf_path: str | Path,
    png_path: str | Path,
    *,
    target_name: str,
    query_name: str,
    window_size: int,
) -> dict[str, int]:
    preflight_plot()
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    panels = _panels(rows)
    page_count = panel_count = 0
    with PdfPages(pdf_path, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        while True:
            page = list(islice(panels, MAX_PANELS_PER_PAGE))
            if not page and page_count:
                break
            count = len(page)
            columns = 1 if count < 4 else 2
            plot_rows = max(1, math.ceil(count / columns))
            figure, axes = plt.subplots(
                plot_rows, columns, figsize=(7 * columns, 3.2 * plot_rows), squeeze=False
            )
            try:
                flat_axes = list(axes.flat)
                if not page:
                    flat_axes[0].text(
                        0.5,
                        0.5,
                        "No target sequence contained enough genes for a complete window.",
                        ha="center",
                        va="center",
                    )
                for axis, (target, series) in zip(flat_axes[:count], page, strict=True):
                    color_map = plt.get_cmap("turbo", max(1, len(series)))
                    for index, (query, x, y) in enumerate(series):
                        axis.plot(x, y, linewidth=1.5, label=query, color=color_map(index))
                    axis.set_title(f"{target_name}: {target}", loc="left", fontweight="bold")
                    axis.set_xlabel(f"Target gene rank (window={window_size})")
                    axis.set_ylabel("Retention (%)")
                    axis.set_ylim(-2, 102)
                    axis.grid(axis="y", alpha=0.25)
                    axis.legend(title=query_name, fontsize="small", frameon=False, ncol=2)
                for axis in flat_axes[count:]:
                    axis.axis("off")
                page_count += 1
                panel_count += count
                figure.suptitle(
                    f"kfFractBias: {target_name} vs {query_name} (page {page_count})",
                    fontweight="bold",
                )
                figure.tight_layout()
                pdf.savefig(figure, bbox_inches="tight")
                if page_count == 1:
                    figure.savefig(png_path, dpi=180, bbox_inches="tight")
            finally:
                plt.close(figure)
            if not page:
                break
    return {
        "pdf_pages": page_count,
        "png_page": 1,
        "panels": panel_count,
        "max_panels_per_page": MAX_PANELS_PER_PAGE,
        "max_queries_per_panel": MAX_QUERIES_PER_PANEL,
    }

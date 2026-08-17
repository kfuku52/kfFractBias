from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path
from typing import Any


def plot_windows(
    rows: list[dict[str, Any]],
    pdf_path: str | Path,
    png_path: str | Path,
    *,
    target_name: str,
    query_name: str,
    window_size: int,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("Plotting requires matplotlib; rerun with --no-plot to omit figures") from exc

    grouped: dict[str, dict[str, list[dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        grouped[str(row["target_seqid"])][str(row["query_seqid"])].append(row)

    if not grouped:
        figure, axis = plt.subplots(figsize=(9, 3))
        axis.axis("off")
        axis.text(
            0.5,
            0.5,
            "No target sequence contained enough genes for a complete window.",
            ha="center",
            va="center",
        )
    else:
        count = len(grouped)
        columns = 1 if count < 4 else 2
        plot_rows = math.ceil(count / columns)
        figure, axes = plt.subplots(plot_rows, columns, figsize=(7 * columns, 3.2 * plot_rows), squeeze=False)
        flat_axes = list(axes.flat)
        for axis, target_seqid in zip(flat_axes, sorted(grouped)):
            for query_seqid, query_rows in sorted(grouped[target_seqid].items()):
                query_rows.sort(key=lambda row: int(row["start_rank"]))
                x_values = [
                    (int(row["start_rank"]) + int(row["end_rank"])) / 2 for row in query_rows
                ]
                y_values = [float(row["retention_percent"]) for row in query_rows]
                axis.plot(x_values, y_values, linewidth=1.5, label=query_seqid)
            axis.set_title(f"{target_name}: {target_seqid}", loc="left", fontweight="bold")
            axis.set_xlabel(f"Target gene rank (window={window_size})")
            axis.set_ylabel("Retention (%)")
            axis.set_ylim(-2, 102)
            axis.grid(axis="y", alpha=0.25)
            axis.legend(title=query_name, fontsize="small", frameon=False, ncol=2)
        for axis in flat_axes[count:]:
            axis.axis("off")

    figure.suptitle(f"kfFractBias: {target_name} vs {query_name}", fontweight="bold")
    figure.tight_layout()
    figure.savefig(pdf_path, bbox_inches="tight")
    figure.savefig(png_path, dpi=180, bbox_inches="tight")
    plt.close(figure)

"""Measure synthetic dense output without checking biological data into Git."""

from __future__ import annotations

import argparse
import json
import resource
import sys
import tempfile
import time
import tracemalloc
from pathlib import Path

from kffractbias.analysis import AnalysisConfig, calculate_fractionation_bias


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", type=int, default=10000)
    parser.add_argument("--queries", type=int, default=100)
    parser.add_argument("--collect-rows", action="store_true")
    parser.add_argument(
        "--trace-memory", action="store_true", help="Measure Python allocations (slows the run)"
    )
    args = parser.parse_args()
    if args.genes < 2 or args.queries < 1:
        parser.error("require at least two genes and one query sequence")
    with tempfile.TemporaryDirectory(prefix="kffractbias-benchmark-") as temporary:
        root = Path(temporary)
        target, query, anchors = (
            root / name for name in ("target.bed", "query.bed", "pairs.anchors")
        )
        target.write_text(
            "".join(
                f"chr{i // (args.genes // 2)}\t{i * 10}\t{i * 10 + 3}\tt{i}\n"
                for i in range(args.genes)
            )
        )
        query.write_text("".join(f"qchr{i}\t0\t3\tq{i}\n" for i in range(args.queries)))
        anchors.write_text("".join(f"t{i}\tq{i % args.queries}\t10\n" for i in range(args.genes)))
        if args.trace_memory:
            tracemalloc.start()
        start = time.perf_counter()
        # This also runs against older revisions for controlled comparisons.
        options = (
            {"collect_rows": args.collect_rows}
            if "collect_rows" in AnalysisConfig.__dataclass_fields__
            else {}
        )
        result = calculate_fractionation_bias(
            AnalysisConfig(
                anchors,
                "jcvi",
                target,
                query,
                root / "out",
                window_size=100,
                make_plot=False,
                **options,
            )
        )
        elapsed = time.perf_counter() - start
        peak = tracemalloc.get_traced_memory()[1] if args.trace_memory else None
        if args.trace_memory:
            tracemalloc.stop()
        summary = json.loads(result.summary_path.read_text())
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if sys.platform != "darwin":
            rss *= 1024
        print(
            json.dumps(
                {
                    "target_genes": args.genes,
                    "query_sequences": args.queries,
                    "gene_rows": summary["counts"]["gene_table_row_count"],
                    "window_rows": summary["counts"]["window_table_row_count"],
                    "elapsed_seconds": round(elapsed, 3),
                    "peak_rss_mib": round(rss / 1024**2, 2),
                    "peak_python_mib": round(peak / 1024**2, 2) if peak is not None else None,
                    "tsv_sha256": {
                        key: summary["output_sha256"][key] for key in ("genes", "windows")
                    },
                },
                indent=2,
                sort_keys=True,
            )
        )


if __name__ == "__main__":
    main()

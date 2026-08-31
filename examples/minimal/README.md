# Minimal precomputed-synteny example

Run from the repository root after installing kfFractBias:

```bash
kffractbias calculate \
  --synteny examples/minimal/pairs.anchors \
  --format jcvi \
  --target-bed examples/minimal/target.bed \
  --query-bed examples/minimal/query.bed \
  --window-size 2 \
  --output-dir example-output \
  --prefix minimal \
  --no-plot
```

The BED files use zero-based, half-open coordinates. The JCVI-style anchors
contain one target and one query identifier per non-comment row.

If you used uv without activating `.venv`, run the command with
`uv run --no-sync kffractbias` instead. This example requires only the base
package; it does not use JCVI, an aligner, or plotting dependencies.

Expected results, excluding the header rows:

- `minimal.genes.tsv`: 8 rows (4 target genes × 2 query sequences).
- `minimal.windows.tsv`: 6 rows (3 complete windows × 2 query sequences).
- `minimal.summary.json`: 3 retained synteny pairs and `schema_version: 3`.
- No PDF/PNG because `--no-plot` is specified.

| Target genes in window | Ranks, inclusive | chrA retained / 2 | chrB retained / 2 |
| --- | --- | --- | --- |
| t1, t2 | 1–2 | 2 / 2 = 100% | 0 / 2 = 0% |
| t2, t3 | 2–3 | 1 / 2 = 50% | 1 / 2 = 50% |
| t3, t4 | 3–4 | 0 / 2 = 0% | 1 / 2 = 50% |

The TSVs are tab-delimited UTF-8 text. See [formats](../../docs/formats.md)
for their columns and [calculation rules](../../docs/methods.md) for the
denominator and window definitions. Omitting `--window-size 2` uses the default
100, which produces no complete windows for this four-gene example.

A successful rerun with the same output directory and prefix replaces the
earlier results. Choose another prefix to retain both runs. For CDS/GFF input,
continue with the [annotation tutorial](../annotations/README.md).

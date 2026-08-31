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

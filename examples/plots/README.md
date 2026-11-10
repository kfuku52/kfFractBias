# Reproduce the README plots

Run from the repository root with an active environment that has the `plot`
extra installed. The [input generator](../../scripts/generate_plot_example.py)
creates two synthetic reference chromosomes and two synthetic queries with two
local copies each. The commands use precomputed anchors.

```bash
work=$(mktemp -d)
python scripts/generate_plot_example.py --output-dir "$work/inputs"

kffractbias calculate \
  --format jcvi \
  --synteny "$work/inputs/query-a.anchors" \
  --target-bed "$work/inputs/reference.bed" \
  --query-bed "$work/inputs/query-a.bed" \
  --target-name Reference --query-name "Query A" --target-seqids chr1 \
  --window-size 40 --step-size 5 \
  --output-dir "$work/results" --prefix retention

printf '%s\n' '{"legend": "graphical", "height_pt": 420, "font_size": 9, "dpi": 180}' \
  > "$work/style.json"
kffractbias report \
  --manifest "$work/inputs/comparisons.tsv" \
  --window-size 40 --step-size 5 --denominator all \
  --plot-config "$work/style.json" \
  --output-dir "$work/results" --prefix shared-reference
```

The [README](../../README.md#example-plots) shows `retention.plot.png` and
`shared-reference.raw.png` from `$work/results`. The latter sums copy detection
in windows that are callable in both queries. The same run also writes
normalised plots and a PDF with profiles for every reference chromosome. See the
[comparison guide](../../docs/comparison.md) for the report settings.

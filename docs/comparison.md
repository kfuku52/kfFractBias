# Shared-reference copy profiles

`report` compares independently assigned local copies on identical reference BED coordinates. It consumes closed JCVI anchors, not sequence alignments or a SynMap JSON export.

```bash
kffractbias report --manifest comparisons.tsv --output-dir report --prefix comparison \
  --window-size 100 --step-size 1 --denominator all
```

The TSV manifest has `analysis_id`, `species`, `expected_copies`, `target_bed`, `query_bed`, `synteny_file`, and `mapping_file`. Paths are relative to the manifest. IDs must be unique. Every target BED must contain identical genes, sequence identifiers, coordinates, strands and order.

A mapping TSV has `gene_id`, `group_id`, `copy`, `segment_id`. Each local group must contain every expected copy label. Assignment must be independent of the measured retention outcome. Unassigned query genes are excluded and counted in the audit. Query sequences and contiguous assignment segments remain separate; no physical joining or genome-wide copy phasing is inferred.

Within each native JCVI block and assignment segment, at least 10 distinct reference anchors define a callable interior, trimmed by five reference-gene ranks at each end. A reference locus is callable only if exactly one local group has evidence for every expected copy. Missing or ambiguous groups remain unknown. Absence of an anchor inside a callable interior is zero detection. Fragment hits are combined by logical OR per local copy before summing.

All species share one callable mask. Complete windows never cross reference sequences. `all` uses every reference gene in window construction; `syntenic` selects genes detected in every analysis before ranking. `callable` selects the common callable loci before constructing windows, regardless of their detected copy count. Thus every callable-based 100-gene window has exactly 100 shared opportunities, but can span more than 100 original reference ranks. Start/end ranks always retain the original reference coordinates. These plots describe the callable subset, not the whole reference genome. The default `all` requires all window loci to be jointly callable. `--minimum-callable-fraction` explicitly permits partial coverage, with only shared callable genes in the denominator. Report this threshold and the coverage distribution when using partial windows. Unknown loci never become zeros.

Outputs are `<prefix>.loci.tsv`, `.report_windows.tsv`, `.report.json`, `.normalised.{pdf,png,svg}`, `.raw.{pdf,png,svg}`, and `.profiles_all.pdf`. Tables include explicit status, shared callable counts, copy values and segment IDs. JSON records source/input/output hashes, settings and coverage audits. `--no-plot` skips figures. Publication uses per-prefix locking, input snapshots and rollback.

When no reference sequence contains enough eligible genes for a complete window, the locus table and summary are still published. The window table contains its header and zero rows, and figures explicitly show "Not estimable". This does not relax callability or window size.

The black trace shows total detected copies divided by expected depth in normalised figures. Coloured traces show individual local-copy detection without dividing by depth. Raw figures retain the summed percentage, which can exceed 100. Boxplots describe overlapping windows and do not constitute independent statistical observations. All plots use common vertical limits across species. Copy traces break at native segment or assignment changes.

`--plot-config` accepts JSON keys `font_family`, `font_size`, `font_files`, `width_pt`, `height_pt`, `dpi`, `reference_seqids`, and `representative_seqids`. Defaults are DejaVu Sans, 8 pt, 518.4 pt wide, 180 pt per species, 300 dpi, and the first three reference sequences. `reference_seqids` explicitly orders displayed reference chromosomes and preserves empty chromosomes as not estimable. It must contain every observed window sequence and may only name sequences present in the reference BED. Font-file paths are relative to this JSON. Text bounds and resolved fonts are verified before publication.

`--raw-unmasked` is a separate descriptive historical mode: it counts native query sequences (or provided copy mappings) without callable-interior masking. It must not be interpreted as a matched absence analysis. Modern JCVI anchors do not reproduce an older CoGe/SynMap analysis merely because window parameters agree.

# Calculation rules

See [formats](formats.md) for input and output fields and the
[minimal example](../examples/minimal/README.md) for values you can reproduce.

**What is counted**

One selected BED row is one counting unit. For CDS/GFF comparison, use one
representative CDS per gene locus. `--isoform-policy error` is the default;
`longest` selects the longest input CDS per known locus, and `all` counts every
mapped identifier. Relationships come from GFF3 Parent or GTF gene_id, not
identifier spelling. The summary reports unresolved loci and the counting unit.

For a target gene and a query sequence, retention is 1 if at least one retained
synteny pair links that gene to any query gene on that sequence; otherwise it is
0. Two query genes on the same sequence still contribute only 1. A target gene
can independently have retention 1 on several query sequences.

`calculate` uses the supplied pairs. It does not perform alignment, chaining,
quota screening, or diagonal filtering. Supply anchors already screened for
your analysis. `compare` and `selfcompare` generate and screen anchors first.
Profiles describe gene retention; the program does not calculate a statistical
significance test or reconstruct missing ancestral loci.

**Sequence selection and denominators**

`--target-seqids` and `--query-seqids` accept comma-separated lists and can be
repeated. `selfcompare --seqids` selects the same set on both axes. Exclusion
regexes apply to both axes. Unknown sequence names or an empty selected genome
are errors. Selection controls profiling; comparison alignment still uses the
prepared genome inputs, so selecting fewer sequences does not itself reduce
the alignment workload.

All synteny data rows are parsed against the complete BED identifier sets
before sequence filtering. A malformed or unknown-ID row outside the selected
sequences is therefore still an error.

- `--denominator all` keeps all selected target genes, including unmatched ones.
- `--denominator syntenic` removes target genes with no retained match to any
  included query sequence, then forms the windows from the remaining genes.

Genes are ordered within each target sequence by BED start, end, then identifier.
Sequence names use deterministic natural ordering. Ranks start at 1 and are
assigned after the denominator filter. If only t1 and t3 of the ordered genes
t1,t2,t3,t4 are syntenic, their syntenic ranks are 1 and 2. Rank 2 in that run
does not refer to the same gene as rank 2 with denominator all.

By default only query sequences represented in retained pairs are output.
`--include-unmatched-query-seqids` includes every selected query sequence.
Explicit `--query-seqids` also includes unmatched members of the specified set,
even without that flag; `selfcompare --seqids` has the same effect. At least
one synteny pair must remain after selection. Selecting only unmatched
sequences gives an error, not an all-zero successful result.

**Windows**

`--window-size W` defaults to 100 genes and `--step-size S` to 1 gene. These are
counts of analyzed BED rows, not nucleotide distances. On each target sequence,
window start ranks are 1, 1+S, 1+2S, and so on. End ranks are start+W-1 and are
inclusive. A window is output only if all W genes exist on that sequence.

No window crosses a sequence boundary. Sequences with fewer than W analyzed
genes contribute no windows, and incomplete trailing windows are omitted.
S greater than W leaves gaps between windows. Gene-level output is still
written when no complete window exists anywhere; that is a successful run.

For each window and query sequence:

```text
retained_count = sum of the W binary per-gene retention values
retention_fraction = retained_count / W
retention_percent = 100 * retention_fraction
```

The number of matches on other query sequences does not change this denominator.
`window_index` restarts at 1 for each target/query sequence combination.
Use start_gene/end_gene to compare intervals across runs with different filters.

**Output size and memory**

Let Q be the number of output query sequences, G the number of analyzed target
genes, and Gc the analyzed gene count on target sequence c. Excluding headers:

```text
gene rows   = G * Q
window rows = Q * sum(max(0, floor((Gc - W) / S) + 1) for each c)
```

TSVs retain zero rows and form dense products over the included sequences.
Streaming controls memory, not row count or disk size. Increasing step size
reduces window rows, but not gene rows. Narrowing sequence selection can reduce
both. The pre-alignment estimate is an upper bound because final retained
pairs and denominator filtering may reduce the actual output.

CLI runs stream rows. The Python API collects result rows by default; use
`AnalysisConfig(collect_rows=False, ...)` for streaming without returned row
tuples. PDF plots include all panels, at most 12 query sequences per panel and
6 panels per page. A target sequence can appear in several panels. PNG is a
preview of the first PDF page only.

**Quota and self retention**

For pairwise comparison, `--quota 1:2` represents one target region to two query
regions; reverse the quota when reversing genomes. This screens block depth,
not an exact per-gene copy-number requirement. Choose the quota from the
biological comparison rather than from the desired retention percentages.

Self-comparison removes identity and mirrored pairs, then counts each retained
nonidentity pair in both directions. Its diagonal bound is a gene-rank distance
within the same chromosome: with the default 300, distances below 300 are
excluded and exactly 300 is allowed. Interchromosomal pairs are not subject to
that distance test. The filter is applied before both chaining and liftover;
upstream identity, C-score, and tandem filtering remain independent.

`selfcompare --depth N` constrains overlapping blocks across both arms of the
same genome. Block intervals use inclusive gene ranks without an overlap
tolerance. It does not allow a region to exceed depth N by appearing on
different axes in different blocks.

The result is `self_synteny_retention`, conditional on extant annotated genes.
It cannot observe an ancestral locus when every descendant copy was lost, so
it has a different interpretation from an outgroup-based pairwise profile.
Reusing self anchors with `calculate --self` requires identical selected BED
gene sets and coordinates on both axes. It removes identities/mirrors and
counts both directions, but does not apply a new depth or diagonal bound.

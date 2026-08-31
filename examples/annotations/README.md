# Synthetic annotation tutorial

Run from the repository root using the [comparison environment](../../README.md)
and LAST on PATH. This tutorial generates small artificial nucleotide sequences
and GFF3 hierarchies; they are test data, not biological genes or a benchmark
for real genomes. The random sequences are not guaranteed to encode proteins.

The default target has 8 loci on target_chr. The query has 8 corresponding loci
on each of query_a and query_b, with sequence divergence between the two copies.
Every locus has one mRNA/CDS representative and an explicit gene Parent chain.
No external data download or CoGe access is needed.

```bash
python scripts/generate_example.py --output-dir annotation-inputs
kffractbias validate \
  --target-cds annotation-inputs/target.cds.fa \
  --target-gff annotation-inputs/target.gff3 \
  --query-cds annotation-inputs/query.cds.fa \
  --query-gff annotation-inputs/query.gff3
```

The generator refuses an existing output directory. Pick a new directory to
generate another copy; it does not overwrite input files. Without environment
activation, prefix the Python and kffractbias commands with `uv run --no-sync`.
Validation should map all 8 target and 16 query identifiers, with no collapsed
isoforms or unresolved loci. Validation itself needs no aligner or JCVI.

**Pairwise comparison**

```bash
kffractbias compare \
  --target-cds annotation-inputs/target.cds.fa \
  --target-gff annotation-inputs/target.gff3 \
  --query-cds annotation-inputs/query.cds.fa \
  --query-gff annotation-inputs/query.gff3 \
  --quota 1:2 --window-size 4 \
  --output-dir annotation-results --prefix pairwise --no-plot
```

With the documented dependencies and default filtering, expect 16 retained
pairs, 16 gene rows, and 10 window rows. Each query sequence retains all four
target genes in every window, giving 100% retention. The selected anchors are
`annotation-results/pairwise.synteny/target.query.lifted.1x2.anchors`.

**Self-comparison of the duplicated query**

```bash
kffractbias selfcompare \
  --cds annotation-inputs/query.cds.fa \
  --gff annotation-inputs/query.gff3 \
  --depth 1 --window-size 4 \
  --output-dir annotation-results --prefix self --no-plot
```

Expect 8 undirected interchromosomal pairs and 16 directed pairs, 32 gene rows,
and 20 window rows. Each query_a window has 100% retention on query_b and 0%
on query_a; query_b has the converse profile. The default diagonal bound of
300 does not remove these interchromosomal pairs, even though both sequences
have only 8 genes. The selected anchors are
`annotation-results/self.synteny/self.self.lifted.1x1.anchors`.

These expectations are sanity checks for the toy inputs, not guarantees of
pair counts on biological data. See [calculation rules](../../docs/methods.md)
for the interpretation of self retention and complete windows.

To exercise BLAST+, install both `blastn` and `makeblastdb`, add `--aligner blast`
to either comparison, and use a different prefix such as pairwise_blast or
self_blast. The default `blastn` task is used; the expected toy pair sets and
counts are the same. Removing `--no-plot` adds PDF/PNG output.

Compare/selfcompare refuse to replace an existing synteny work directory unless
`--force` is given. For a fresh run that preserves prior results, use a new
prefix. To change only the windows, reuse the BED/anchors with `calculate` as
shown in the [README](../../README.md), adding `--self` for self anchors.

# kfFractBias

`kfFractBias` is a Python 3 command-line program for calculating and plotting
gene-retention and fractionation-bias profiles. It is a maintained fork and
offline redesign of
[`SynMapFractBiasAnalysis`](https://github.com/bjoyce3/SynMapFractBiasAnalysis).

The original FractBias implementation was developed by Blake L. Joyce, Asher
Haug-Baltzell, Sean Davey, Matthew Bomhoff, James C. Schnable, and Eric Lyons.
This repository contains the maintained Python 3 implementation; the original
Python 2 and notebook sources remain available in the Git history.

## What is different

- Python 3 package and the `kffractbias` executable
- no CoGe account, genome ID, API, or JWT requirement
- calculation from either JCVI anchors or legacy SynMap/DAGCHAINER output
- end-to-end local comparison from target/query CDS FASTA and GFF annotations
- within-genome self-synteny retention with chromosome-aware JCVI chaining
- local JCVI MCscan and QUOTA-ALIGN execution
- explicit target:query syntenic-depth quota
- long-form gene and sliding-window TSV outputs
- JSON provenance with input hashes and parameters
- non-interactive PDF and PNG plots

After the software and its dependencies are installed, analysis does not need
network access.

## Installation

Python 3.11 or newer on Linux or macOS is required. Output locking and replacement
are supported on local POSIX filesystems; Windows and network filesystems are
not supported or tested.

```bash
python -m pip install '.[plot]'
kffractbias version
```

The base installation has no third-party runtime dependencies and supports
`calculate --no-plot`, `validate`, `formats`, and `version`. Use `.[plot]` for
plots or `.[compare]` for plots plus the end-to-end JCVI commands. For a fully
locked environment, install [uv](https://docs.astral.sh/uv/) and run
`uv sync --locked --extra all`.

The `compare` extra includes JCVI, but its default `last` aligner also requires
the external `lastal` and `lastdb` executables. JCVI QUOTA-ALIGN requires its
mixed-integer solver dependency.
For `--aligner blast`, install NCBI BLAST+ (`blastn` and `makeblastdb`) instead.
The search task is explicitly `blastn` for divergent CDS; choose
`--blast-task dc-megablast` or `--blast-task megablast` if those search heuristics
are appropriate. The latter can miss short or regularly diverged sequences
without a long exact seed. The task and all external commands are recorded;
nonzero BLAST exit codes stop the run before JCVI processing.
Executable, JCVI/solver, plot, mapping, sequence-selection, and output-path
checks run before alignment. A dense output row-count upper bound is printed
and saved as `preflight.json` in the synteny directory.

## Commands

### Calculate from precomputed synteny

`calculate` accepts a JCVI `.anchors` file or SynMap/DAGCHAINER output. The
target and query BED identifiers must match the synteny identifiers. BED
coordinates use the standard zero-based, half-open convention. Target and
query gene identifier sets must be disjoint for pairwise analysis so pair
orientation cannot be guessed incorrectly.

```bash
kffractbias calculate \
  --synteny target.query.lifted.1x2.anchors \
  --format jcvi \
  --target-bed target.bed \
  --query-bed query.bed \
  --target-name Sorghum_bicolor \
  --query-name Zea_mays \
  --window-size 100 \
  --denominator all \
  --output-dir results \
  --prefix sorghum_maize
```

### Compare two annotations end-to-end

`compare` detects matching GFF feature/attribute identifiers, prepares JCVI
BED and CDS inputs, runs MCscan and QUOTA-ALIGN, and then calculates the
fractionation profile.

```bash
kffractbias compare \
  --target-cds target.cds.fa \
  --target-gff target.gff3 \
  --query-cds query.cds.fa \
  --query-gff query.gff3 \
  --target-name Sorghum_bicolor \
  --query-name Zea_mays \
  --quota 1:2 \
  --window-size 100 \
  --cpus 8 \
  --output-dir results \
  --prefix sorghum_maize
```

The quota is always explicit: `1:2` means one expected target region for two
expected query regions. Reversing target and query also requires reversing the
quota.

CDS-to-GFF mapping is strict by default: every FASTA identifier must map to the
selected GFF feature and attribute. If incomplete mapping is intentional, set
`--minimum-mapping-fraction` explicitly; the matched and total counts are
recorded in the summary.

Use one representative CDS per gene locus. When `Parent` in GFF3 or `gene_id`
in GTF identifies multiple mapped isoforms of a locus, the default
`--isoform-policy error` rejects the input. Choose `--isoform-policy longest`
to select the longest input CDS per locus; equal lengths use deterministic
identifier order. `--isoform-policy all` explicitly counts each mapped
identifier, including separate transcripts of the same locus. The mapping
fraction is checked **before** isoform selection, so intentional selection does
not look like missing annotation.

The selected identifiers remain unchanged in BED and anchors. Windows and
denominators count these selected BED rows. Parent relationships are resolved
through transcript/CDS hierarchies, not inferred from identifier spelling.
Where the annotation supplies no gene relationship, an identifier remains a
separate unit; the summary reports `unresolved_locus_count`. In that case,
provide representative CDS inputs yourself to ensure a gene-locus analysis.
The summary also records the policy, original match count, selected count,
collapsed isoform count, and counting unit. Representative-only inputs retain
their previous identifiers and counts.

### Compare a genome to itself

`selfcompare` measures within-genome retention across non-redundant self-synteny
blocks. It uses JCVI for alignment, C-score/tandem filtering, chaining, and
liftover. A local adapter removes identity and mirrored hits and applies the
diagonal distance **only within the same chromosome**, before both chaining
and liftover. Interchromosomal pairs are independent of concatenated BED ranks
and unrelated chromosomes.

```bash
kffractbias selfcompare \
  --cds polyploid.cds.fa \
  --gff polyploid.gff3 \
  --name Polyploid_species \
  --depth 1 \
  --window-size 100 \
  --cpus 8 \
  --output-dir results \
  --prefix polyploid_self
```

`--depth N` applies a symmetric `N:N` maximum block depth. Overlap constraints
combine both arms on the same genome, including regions that appear on the
query side of one block and the target side of another. Block intervals use
inclusive gene ranks without an overlap tolerance. It must be chosen
from the expected homeologous structure rather than inferred automatically.
The default `--self-hit-percent 98` is JCVI's near-self identity cutoff and can
be changed for unusually recent polyploidy. The default `--diagonal-bound 300`
excludes intrachromosomal anchors whose gene-rank difference is **less than**
300; a difference of exactly 300 is allowed. Change that option for a different
intrachromosomal diagonal definition. There is no additional hidden 40- or
300-rank cutoff in the adapter or quota screen. JCVI's upstream alignment and
tandem filtering still apply independently.

The result is explicitly recorded as `self_synteny_retention`. It measures
retention asymmetry conditional on genes present in the annotated genome; it
cannot observe ancestral loci for which every descendant copy was lost and is
therefore not equivalent to outgroup-based fractionation bias from `compare`.

### Validate annotations

```bash
kffractbias validate \
  --target-cds target.cds.fa --target-gff target.gff3 \
  --query-cds query.cds.fa --query-gff query.gff3
```

Use `--target-feature` with `--target-attribute`, and the corresponding query
options, when automatic GFF identifier detection is not appropriate.
`validate` also checks FASTA headers, duplicate identifiers, nonempty sequences,
and nucleotide IUPAC characters. Wrapped, lowercase, and gzip-compressed FASTA
are accepted; it does not validate translation, reading frame, or annotation
biological correctness. It supports the same `--isoform-policy` as comparison.

Every non-comment synteny row is validated. A malformed row or a row that does
not contain one target and one query identifier stops the analysis with its
line number instead of being silently discarded. Duplicate valid pairs are
deduplicated and counted in the summary.

Other commands are `kffractbias formats` and `kffractbias version`.

## Outputs

For a prefix such as `sorghum_maize`, kfFractBias writes:

- `sorghum_maize.genes.tsv`: retention and matching query IDs for each target
  gene and query sequence combination
- `sorghum_maize.windows.tsv`: sliding-window retention fractions and percent
- `sorghum_maize.summary.json`: parameters, counts, input hashes, and outputs
- `sorghum_maize.plot.pdf`: all chromosome-wise profiles, paginated at six
  panels per page and twelve query sequences per panel
- `sorghum_maize.plot.png`: a preview of the first PDF page; the summary records
  page/panel counts (use the PDF for all profiles)
- `sorghum_maize.synteny/`: retained JCVI working data from `compare`

`selfcompare` writes the same output set and records identity/mirror filtering,
interchromosomal and intrachromosomal pair counts, symmetric depth, and the
interpretation limitation in the JSON summary.

`--denominator all` uses every target gene in each window. `--denominator
syntenic` first removes target genes without any retained query match, matching
the two denominator choices exposed by the original FractBias implementation.

By default, output tables include only query sequences represented by retained
synteny pairs. This avoids a target-gene by every-query-contig cross product on
fragmented assemblies. Use `--include-unmatched-query-seqids` when explicit
zero-retention profiles for every selected query sequence are required.

The summary schema records input and output SHA-256 hashes, input and selected
gene counts, synteny record and duplicate counts, Python/platform/package
versions, and JCVI/aligner versions when available. Inputs are hashed before
any preparation and verified again before output commit. Each run reads fixed,
private copies of its original inputs; source and prepared BED/CDS hashes are
recorded separately. Dense TSV rows are streamed without removing zero rows.
The library still returns collected rows by default; use
`AnalysisConfig(collect_rows=False, ...)` for CLI-like memory usage.

One per-prefix lock protects preparation, alignment, calculation, and commit.
`--force` permits replacement of existing synteny work only after the new run
succeeds; it never bypasses the lock. Results and synteny work are staged in a
private directory and replaced together with rollback on ordinary failures.
Keep enough disk space for the input copies, new results, and old results.
Locks refuse symlinks and hardlinks and remain on disk after release. They do
not promise power-loss/SIGKILL atomicity across multiple output files.
Individual renames are atomic; readers needing a consistent output set should
wait for completion. The summary is installed last and includes TSV/plot hashes.

For diagnostics, add `--keep-failed-work`. On failure, the CLI prints its unique
`.PREFIX.staging-*` directory containing fixed inputs, partial work, and
`failure.json`. Successful comparisons retain per-command stdout/stderr logs,
commands, exit codes, and elapsed times in `PREFIX.synteny/logs/`; stage timings
are recorded in the summary. Recovery files are always retained if rollback
itself fails. These directories can contain genomic data; remove them manually
after diagnosis when they are no longer needed.

To change the window, sequence selection, or plotting without rerunning an
expensive pairwise alignment, reuse the retained inputs with a new prefix:

```bash
kffractbias calculate \
  --synteny results/sorghum_maize.synteny/target.query.lifted.1x2.anchors \
  --target-bed results/sorghum_maize.synteny/target.bed \
  --query-bed results/sorghum_maize.synteny/query.bed \
  --window-size 50 --step-size 5 \
  --output-dir results --prefix sorghum_maize_window50
```

For existing **self** anchors, add `--self` and supply the same BED on both axes.
This preserves self-retention interpretation, identity/mirror removal, and
bidirectional counting without another alignment. It does not apply a new quota
or diagonal bound; regenerate anchors with `selfcompare` to change those.

```bash
kffractbias calculate --self \
  --synteny results/polyploid_self.synteny/self.self.lifted.1x1.anchors \
  --target-bed results/polyploid_self.synteny/self.bed \
  --query-bed results/polyploid_self.synteny/self.bed \
  --window-size 50 --output-dir results --prefix polyploid_self_window50
```

Without `--self`, `calculate` requires disjoint target/query gene IDs.

A small current-format example is available in [`examples/minimal`](examples/minimal).

## Development

```bash
uv sync --locked --extra test --extra plot
uv run --no-sync python scripts/check.py
```

The opt-in pairwise and self-synteny JCVI/LAST integration tests (plus BLAST+ when installed) run in CI and can be run locally with
`uv sync --locked --extra test --extra all`, then
`uv run --no-sync python scripts/check.py --full --integration`.
The full check builds and installs the wheel without runtime dependencies,
then extracts and tests the sdist in a separate environment. CI caches uv
downloads by Python version and lockfile and preserves required check names.
Use `uv run --no-sync python scripts/benchmark.py --genes 10000 --queries 100`
for a repeatable synthetic workload; `--trace-memory` measures allocations
separately from uninstrumented timings.
See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the full checklist.

## Citation

If you use kfFractBias, cite the original FractBias publication:

> Joyce BL, Haug-Baltzell A, Davey S, Bomhoff M, Schnable JC, Lyons E.
> FractBias: a graphical tool for assessing fractionation bias following
> polyploidy. Bioinformatics. 2017;33(4):552–554.
> <https://doi.org/10.1093/bioinformatics/btw666>

## License

MIT. The original copyright and license are retained in [`LICENSE`](LICENSE).

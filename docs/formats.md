# Input and output formats

Text inputs and TSV outputs use UTF-8. Text input files ending in `.gz` are
read as gzip-compressed text. Blank lines and lines beginning with `#` are
ignored in BED and synteny input. See [calculation rules](methods.md) for the
meaning of retention and the [minimal example](../examples/minimal/README.md)
for a complete runnable dataset.

**CDS FASTA and annotations**

FASTA identifiers are the first whitespace-delimited token after `>`. They
must be unique and have nonempty nucleotide sequences. Wrapped sequences,
lowercase letters, and IUPAC symbols `ACGTURYSWKMBDHVN` are accepted. Empty
headers/records, sequence before the first header, and other characters are
errors. The parser checks syntax, not reading frames, translations, or
biological correctness.

GFF3/GTF uses standard annotation columns with one-based inclusive coordinates;
prepared BED converts them to zero-based half-open coordinates. FASTA IDs must
map to a selected feature/attribute pair. Automatic detection considers mRNA,
transcript, gene, and CDS; explicit feature/attribute overrides must be paired.
Pairwise overrides are `--target-feature`/`--target-attribute` and their query
equivalents; selfcompare uses `--feature`/`--attribute`.

The default mapping fraction is 1.0. The fraction is checked before isoform
selection. GFF3 Parent/GTF gene_id relationships determine known loci; missing
relationships are reported rather than guessed. See the README for the
error/longest/all isoform policies and the
[annotation tutorial](../examples/annotations/README.md) for a synthetic example.

**BED**

BED input is tab-delimited, with at least four columns:

| Position | Field | Requirement |
| --- | --- | --- |
| 1 | sequence ID | Nonempty chromosome/contig name. |
| 2 | start | Integer, zero-based, at least zero. |
| 3 | end | Integer, exclusive endpoint, greater than start. |
| 4 | gene ID | Nonempty, unique within that BED. |
| 5 | score | Optional; not used in retention calculations. |
| 6 | strand | Optional; stored as `.` if absent. |

Additional columns are not used. Multiple exons must not be supplied as
separate BED rows with the same ID. The program orders rows; the input does not
need to be sorted. Pairwise target/query ID sets must be disjoint across the
complete BEDs. With `calculate --self`, selected genes and their sequence,
start, end, and strand must be identical on both axes.

**JCVI anchors**

The first two tab-delimited columns are gene IDs; additional score columns are
not used by `calculate`. Both target-query and query-target orientations are
accepted when membership is unambiguous. Block comments such as `###` do not
affect profiling. Duplicate pairs are counted in the summary and contribute
only once to retention.

```text
###
t1	q1	100
t2	q2	90
```

Every data row must resolve to the complete BED identifier sets, including
rows on sequences that will later be excluded. The parser does not perform
alignment quality or quota screening on supplied anchors.

**Legacy SynMap extended DAGCHAINER**

`--format synmap` accepts SynMap's genomic-coordinate records with embedded
`||` subfields. It does not accept every generic DAGCHAINER layout. After
expanding `||` to tabs, at least 20 fields are required. The one-based candidate
ID positions are 5 and 8 on the left, 17 and 20 on the right. These correspond
to the annotation name or CoGe numeric feature ID in the supported layout.
Each side must resolve unambiguously against the supplied BEDs.

For example, this is one record with literal tabs between the outer fields:

```text
a1_chr1	chr1||1||2||t1||1||CDS||101||1||99	1	2	b2_chrA	chrA||3||4||q1||1||CDS||202||1||99	3	4	1e-20	10
```

It can resolve to BED IDs t1/q1 or 101/202, depending on the BEDs. If both
candidate IDs match distinct genes in the same BED, the row is ambiguous and
is rejected. The names and numeric IDs are not obtained from a CoGe API.

`--format auto` inspects the first data row: `||` selects SynMap, otherwise
JCVI. Use an explicit format for known input. `kffractbias formats` prints a
short format summary.

**genes.tsv**

Each analyzed target gene has one row per output query sequence, including
zero-retention rows. Column order is fixed:

| Column | Meaning |
| --- | --- |
| target_seqid | Target chromosome/contig. |
| target_gene | Unchanged target BED identifier. |
| target_rank | One-based rank among analyzed genes on this target sequence. |
| query_seqid | Query chromosome/contig. |
| retained | Integer 0 or 1: at least one retained pair on this query sequence. |
| query_genes | Naturally sorted matching query IDs, separated by `;`; empty when retained is 0. |

Rows are ordered by target sequence, target rank, then query sequence. Ranks
refer to the genes remaining after the selected denominator filter. These are
not BED coordinates or stable ranks across differently filtered runs.

**windows.tsv**

| Column | Meaning |
| --- | --- |
| target_seqid | Target chromosome/contig. |
| query_seqid | Query chromosome/contig. |
| window_index | One-based index within this target/query sequence combination. |
| start_rank | One-based rank of the first analyzed gene in the window. |
| end_rank | One-based inclusive rank of the final analyzed gene. |
| start_gene | Identifier at start_rank. |
| end_gene | Identifier at end_rank. |
| retained_count | Sum of binary per-gene retention in this window and query sequence. |
| window_size | Number of analyzed genes in the complete window. |
| retention_fraction | retained_count / window_size, in [0, 1]. |
| retention_percent | 100 × retention_fraction, in [0, 100]. |

Rows are ordered by target sequence, query sequence, then window index.
Fractions and percentages are decimal text with up to ten significant digits;
counts permit independent recalculation. There are no partial-window rows.
A successful analysis can write only the header if no complete windows exist.

**summary.json, schema version 3**

Read `schema_version` before consuming the remaining keys. Package version and
JSON schema version are separate; a development package update need not change
the schema. Consumers should use field names rather than JSON member order.

| Top-level key | Contents |
| --- | --- |
| schema_version | Integer 3 for this schema. |
| program, program_version | Program name and installed package version. |
| analysis_mode | `pairwise_fractionation_bias` or `self_synteny_retention`. |
| target_name, query_name | User-facing genome labels. |
| runtime | Python/platform and installed package versions; unavailable optional packages are null. |
| parameters | Window/step, denominator, resolved and requested synteny formats, resolved sequence lists, exclusion regex, unmatched-query flag. |
| counts | Input, selected, analyzed, pair, and output-row counts described below. |
| inputs | Input labels mapped to original/final prepared paths and raw-file SHA-256 hashes. Source FASTA/GFF and prepared BED/CDS are recorded separately for comparisons. |
| outputs | genes/windows TSV and PDF/PNG paths; plot paths are null when plotting is off. |
| output_sha256 | SHA-256 of the TSV/plot files; plot hashes are null when plotting is off. The summary does not hash itself. |
| timings_seconds | Measured stage times. Nested synteny/preflight stages overlap; do not sum them as elapsed wall time. The successful summary is written before commit, so it does not include commit time. |
| plot | Page/panel counts and limits, including png_page 1; empty when plotting is off. |
| metadata | Comparison-generation details, tools, quota, annotation mapping, and self interpretation where applicable; empty for a basic pairwise calculate run. |

Common counts have these meanings:

| Key in counts | Meaning |
| --- | --- |
| input_target_gene_count, input_query_gene_count | BED rows before sequence selection. |
| target_gene_count, query_gene_count | BED rows after sequence selection, before the target denominator filter. |
| analyzed_target_gene_count | Target rows after the denominator filter. |
| analyzed_query_sequence_count | Number of query sequences actually emitted. |
| input_synteny_record_count | Non-comment/nonblank rows read, including duplicate records. |
| duplicate_synteny_pair_count | Repeated valid pairs collapsed by the parser. |
| sequence_filtered_synteny_pair_count | Unique parsed pairs discarded by sequence selection. |
| synteny_pair_count | Selected pairs; in self mode, unique undirected nonidentity pairs. |
| gene_table_row_count, window_table_row_count | Data rows written, excluding headers. |

Self mode adds input_synteny_pair_count (after sequence selection, before
identity/mirror removal), removed_identity_pair_count,
removed_mirrored_pair_count, directed_synteny_pair_count, and the
intrachromosomal_pair_count/interchromosomal_pair_count split of undirected
retained pairs. Exact duplicate rows are already accounted for by the parser.

Comparison mapping metadata includes feature, attribute, fasta_gene_count,
matched_gene_count (before isoform selection), selected_gene_count,
collapsed_isoform_count, unresolved_locus_count, isoform_policy, and counting_unit.
An unresolved locus means the annotation did not establish a known gene locus;
it is not evidence that the identifier is biologically an independent gene.

**Working files and publication**

Comparisons retain `PREFIX.synteny/`, including prepared BED/CDS, anchors,
`preflight.json` output-row upper bounds, and `logs/` command output. Each
command log has a JSON record with its command, working directory, return code,
elapsed seconds, and log path. Failed work is kept when requested; see
[troubleshooting](troubleshooting.md).

Results are staged and the summary is installed last. Readers should wait for
run completion and check hashes when they need a consistent set; multi-file
atomicity across power loss or SIGKILL is not promised. A successful rerun
replaces existing results at the same prefix, and `--no-plot` removes old plots
for that prefix. Lock-file presence alone does not indicate an active run.

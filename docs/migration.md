# Updating an older checkout

The current `master` line is a development build identified by a `.dev` version.
Historical 0.1.3 and 0.1.4 were untagged milestones; different commits could
report the same version while implementing different behavior. The
[changelog](../CHANGELOG.md) describes those milestones without presenting them
as published releases. Tagged releases, when published, use `vX.Y.Z`.

Keep the old input files, output files, summary JSON, command, and dependency
environment before changing an analysis. Record the exact old commit as well
as the displayed package version. To record the source revision of a checkout:

```bash
git rev-parse HEAD
kffractbias version
```

For uv without activation, use `uv run --no-sync kffractbias version`.
Save the matching `uv.lock` and external aligner versions. Summary JSON records
package/tool versions and input hashes, but does not automatically capture the
Git commit of an editable installation.

| Older behavior or assumption | Current behavior and migration |
| --- | --- |
| Python 3.9+ in the old metadata | Use Python 3.11+ on Linux/macOS with a local POSIX output filesystem. |
| JCVI and matplotlib installed as base dependencies | Select plot or compare/all extras as described in the [README](../README.md). Base calculate needs `--no-plot`. |
| The same gene IDs used in two genomes | Pairwise IDs must be disjoint. Prefix IDs consistently in each genome's FASTA/GFF/BED/anchors; do not rename only one file. |
| Partial CDS-to-annotation mapping accepted | Mapping defaults to 100%. Inspect the missing IDs first; use an explicit minimum fraction only if loss is intentional. |
| Multiple transcript IDs counted as independent genes | Choose representative input or `--isoform-policy longest`; use `all` only when identifier-level counting is intended. |
| Malformed synteny rows silently skipped | Every data row must have valid structure and resolve to the supplied BEDs, even if its sequence is later excluded. |
| Native JCVI self filtering and quota behavior | The chromosome-aware adapter fixes interchromosomal filtering and shares self depth constraints across both axes. Regenerate old self anchors to obtain these corrections. |
| Native BLAST search-task defaults | BLAST searches explicitly use `blastn`; overrides are recorded. Regenerate anchors if the old task missed divergent matches. |
| One PNG represents the entire analysis | PNG shows only the first PDF page. Use the PDF for every chromosome/query panel. |
| Legacy Python 2 scripts, notebooks, and bundled biological datasets in the checkout | They remain in [the pre-removal Git tree](https://github.com/kfuku52/kfFractBias/tree/af8130640d3d097b8038b291b9dd40e6a6d5050e). Current examples use small synthetic data. |

Schema version 3 is described in [formats](formats.md). Downstream programs
should inspect `schema_version`, use column/key names, and distinguish input,
selected, and analyzed counts. Do not infer genomic coordinates from ranks or
assume ranks are stable when changing the denominator or sequence selection.

Use a new result prefix while comparing old and new runs. Compare parameters,
counting units, selected sequence sets, and anchor provenance before attributing
retention changes to biology. Different isoform choices, anchors, or denominator
filters can intentionally change results.

Retained pairwise anchors can be reused with `calculate`; retained self anchors
require `calculate --self` and the same BED on both axes. That reuse changes
profiling only. It cannot repair anchors lost or incorrectly selected by an
older aligner, self filter, or quota stage. See the reuse commands in the
[README](../README.md) and [calculation rules](methods.md).

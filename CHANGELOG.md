# Changelog

All notable changes are documented here. This project follows Semantic
Versioning after the first tagged release. Earlier 0.1.x headings identify
untagged development milestones, not published release artifacts. See the
[migration guide](docs/migration.md) when updating an older checkout.

## Unreleased

- Preserve literal GTF attribute values and quoted notes, reject malformed
  attributes with source lines, and prevent gene/transcript name collisions
  from merging independent loci during isoform selection.
- Count overlapping arms of the same self-synteny block separately in depth
  constraints, with interval-coverage oracle and default-bound regressions.
- Reject gene IDs that collide with comments or ID-list/format delimiters,
  with consistent FASTA, mapped GFF, and BED diagnostics.
- Add `validate --pairwise`, exact output row counts and `--max-output-rows`,
  and successful CLI elapsed time including commit and cleanup.
- Adapt JCVI 1.6.6 quota help strings for Python 3.14 argparse and exercise
  real LAST/BLAST comparisons on Python 3.12 and 3.14 in CI.

- Reject empty/invalid GFF and BED sequence IDs and strand values, and reject
  opposite-strand segments for one mapped identifier with file/line diagnostics.
  Preserve valid multi-segment CDS merging and normalize unknown strands.
- Record import-time Git revision/dirty state and a Python-source fingerprint
  in `runtime.source`; preserve source identity in Git-free wheel/sdist installs.
- Check maintained Markdown examples and internal links automatically, run
  tutorial examples in the existing numeric/integration checks, and add macOS
  coverage plus comparison-extra installation on Python 3.14.
- Identify the source development line as `0.2.0.dev0`, distinct from the old
  untagged 0.1.4 implementation; keep citation and lockfile versions in sync.
- Complete source-install and release instructions, repair dependency guidance,
  and document window/rank semantics, output replacement, formats, and migration.
- Add reproducible annotation examples and expected retention values, include
  documentation in the sdist, and enable private vulnerability reporting.
- Fix interchromosomal self-synteny loss at the default diagonal bound and
  enforce self quota across both arms of the same genome.
- Specify the BLAST nucleotide search task explicitly, support task overrides,
  check BLAST exit codes, and test known pair sets with both LAST and BLAST+.
- Hold the output lock across the entire comparison; use private input
  snapshots and work directories, safe lock files, and recoverable replacement
  of work and results together.
- Validate malformed/empty nucleotide FASTA and add explicit error/longest/all
  isoform policies with counting-unit metadata.
- Stabilize natural-sort ties across hash seeds; stream dense tables and
  paginate plots without discarding output rows or PDF panels.
- Add early dependency/selection checks, retained failure diagnostics and
  command timings, isolated distribution checks, and a complete sdist manifest.
- Reuse self-synteny anchors with `calculate --self` when changing windows,
  selected sequences, or plotting without repeating alignment.
- Cache locked CI dependencies and add fast development checks, independent
  numeric oracles, scientific integration assertions, and a synthetic benchmark.

## 0.1.4 - 2026-08-19 (untagged development milestone)

- Reject input/output path collisions and unsafe `--force` work-directory use.
- Validate every synteny row and reject ambiguous pairwise gene identifiers.
- Stage output sets transactionally, lock concurrent prefixes, and record
  input/output hashes plus runtime versions.
- Make annotation mapping strict by default and reduce fragmented-assembly
  output expansion.
- Split optional runtime dependencies, add a lock file, and expand CI coverage.
- Add current examples and repository maintenance documentation.

## 0.1.3 (untagged development milestone)

- Provide the maintained offline Python 3 CLI with pairwise and self-synteny
  analysis.

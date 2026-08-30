# Changelog

All notable changes are documented here. This project follows Semantic
Versioning after the first tagged release.

## Unreleased

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

## 0.1.4 - 2026-08-19

- Reject input/output path collisions and unsafe `--force` work-directory use.
- Validate every synteny row and reject ambiguous pairwise gene identifiers.
- Stage output sets transactionally, lock concurrent prefixes, and record
  input/output hashes plus runtime versions.
- Make annotation mapping strict by default and reduce fragmented-assembly
  output expansion.
- Split optional runtime dependencies, add a lock file, and expand CI coverage.
- Add current examples and repository maintenance documentation.

## 0.1.3

- Provide the maintained offline Python 3 CLI with pairwise and self-synteny
  analysis.

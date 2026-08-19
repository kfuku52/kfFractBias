# Changelog

All notable changes are documented here. This project follows Semantic
Versioning after the first tagged release.

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

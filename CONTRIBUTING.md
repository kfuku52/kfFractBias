# Contributing

Clone the repository as described in the [README](README.md), enter its root,
and use Python 3.11 or newer with the locked development environment:

```bash
uv sync --locked --extra test --extra plot
```

Before opening a pull request or pushing directly to `master`, run:

```bash
uv run --no-sync python scripts/check.py
```

For coverage, wheel/sdist build, and isolated distribution tests, use
`uv run --no-sync python scripts/check.py --full`. For changes to synteny,
install the locked `test` and `all` extras and LAST, then add `--integration`.
Install BLAST+ as well to exercise both aligners as CI does. BLAST tests skip
when `blastn` is absent; a local run with skips is not a check of that aligner.
Use `--no-sync` while checking so an invocation does not remove optional
packages from an already prepared integration environment.
When repeating `uv sync` in that environment, retain `--extra all`: selecting
only `plot` removes the comparison dependencies. Initial sync may download
packages; `--offline` can be added when the required packages are already cached.

## Choosing checks

All commands below run from the repository root in the environment above.
For a focused check, use `uv run --no-sync python -m pytest -q -rs` followed
by the test paths in the table. These are starting points during editing;
the standard `scripts/check.py` check is still required before push.

| Change | Focused tests | Additional verification |
| --- | --- | --- |
| CLI options, README, docs, examples | `tests/test_cli.py tests/test_documentation.py` | Real comparison examples require the integration tier below. |
| BED/FASTA/GFF/synteny parsing or isoform selection | `tests/test_io.py tests/test_analysis.py tests/test_preflight.py` | Include malformed, ambiguous, and duplicate input cases. |
| Retention, windows, sequence selection, denominator or streaming | `tests/test_profiles.py tests/test_analysis.py` | Check known numerical results, not just successful output. |
| Output safety, locks, snapshots, provenance | `tests/test_run.py tests/test_safety.py tests/test_provenance.py` | Filesystem behavior also runs on Linux and macOS in CI. |
| Plotting | `tests/test_plotting.py` | Requires the plot extra; skipped tests do not verify figures. |
| JCVI orchestration, self filtering/depth, compatibility adapter | `tests/test_jcvi.py tests/test_selfscan.py tests/test_preflight.py` | Use all extras for JCVI/solver tests, then real integration. |
| Packaging, version, lockfile, distribution contents | `tests/test_metadata.py` | Run `scripts/check.py --full`; follow the version procedure below when applicable. |

Verification tiers:

- **Fast, local:** `uv run --no-sync python scripts/check.py` checks Ruff lint
  and formatting, mypy, and pytest. After environment setup it needs no network,
  real genomes, or external aligners, provided `KFFRACTBIAS_RUN_INTEGRATION`
  is unset. Pytest examples use temporary directories and assert exact results.
  To see skip reasons, use `uv run --no-sync python -m pytest -q -rs`.
- **Distributions:** `uv run --no-sync python scripts/check.py --full` adds
  coverage (80% threshold), wheel/sdist builds, and isolated installation/tests.
  It writes `.coverage` and `dist/`; isolated environments use temporary
  directories. Build/setup can access package indexes unless the caches and
  offline configuration suffice. This is not a real-data analysis.
- **Real alignment:** after `uv sync --locked --extra test --extra all` and
  installing LAST/BLAST+ as described in [troubleshooting](docs/troubleshooting.md),
  run `uv run --no-sync python scripts/check.py --integration` (add `--full`
  when distribution checks are also needed). It uses small synthetic inputs
  in temporary directories, without genomic downloads. Check both LAST
  executables and both BLAST+ executables first. BLAST skips are not success
  for that aligner; missing required tools must be reported, not bypassed.

For a minimal execution check, run
`uv run --no-sync python -m pytest -q -rs tests/test_documentation.py`.
It checks local links and CLI syntax and executes the minimal calculate and
annotation-validation tutorials in temporary directories. Success means pytest
exits zero with those execution tests passing, including the documented
retention values and annotation counts. It does not run every shell fence,
fetch remote links, or exercise real comparison/plot examples.

Benchmarks, real-genome runs, and the network-dependent CI dependency audit
are separate from the fast suite. Do not use them as an automatic smoke check.

## Test and contribution conventions

Changes to synteny generation should also pass the real JCVI/LAST integration
test documented in the README. New parsing behavior must include malformed,
ambiguous, and duplicate-input tests. Never commit biological input data,
generated result directories, credentials, or local environments.

Synteny tests assert known pair sets, depth, diagonal boundaries, and retention
values, not just nonempty files. The fast suite includes explicit retention
expectations for windowing and denominator selection, README CLI parsing,
subprocess hash-seed checks, multiprocess locking, input mutation, source
provenance, documentation examples/links, and
commit/rollback failure injection. The optional JCVI core
tests do not require LAST. Use synthetic fixture generators for new tests.

The fast suite runs the minimal example and annotation validation directly from
their documented commands. Keep runnable commands in bash/sh/shell/console
fences so they can be checked. The annotation integration tests use the same
generator and documented comparison commands, including the optional plot and
BLAST variants, and check known pairs, table sizes, and retention values.

CI preserves the existing Linux/Python 3.11–3.14 check names and adds a
macOS/Python 3.12 fast job for filesystem/locking behavior. The Linux/Python 3.14
job installs the comparison extras as well; external LAST/BLAST integration
runs on Linux/Python 3.12 and 3.14. Distribution checks compare source fingerprints
between checkout, wheel, and sdist, including a wheel without Git provenance.

Run `python scripts/benchmark.py --genes 10000 --queries 100` to measure dense
streaming tables. Compare equal Python/dependency environments and TSV hashes;
measure Python allocation peaks separately with `--trace-memory`, since the
tracer changes timing. Do not make noisy wall-clock thresholds a CI gate.

The quality job caches the same locked optional environment as integration.
Inspect setup/cache, aligner installation, and test-step durations separately
when assessing CI latency. A cache hit is not a guarantee of a fixed speedup.

Keep `pyproject.toml`, `CITATION.cff`, and `CHANGELOG.md` synchronized when
changing the version. `uv.lock` and installed package metadata must also match;
the metadata test checks consistency rather than a fixed version string.

For a version change or release preparation:

1. Set the intended version in `pyproject.toml` and `CITATION.cff`. Development
   snapshots use a `.devN` suffix. For a release, move the relevant Unreleased
   changes into a dated version section and update the README release-status
   text. Historical untagged milestones must remain labeled as such.
2. Regenerate the lockfile and reinstall this checkout, without upgrading
   unrelated dependencies:

   ```bash
   uv lock
   uv sync --locked --extra test --extra all
   uv run --no-sync kffractbias version
   ```

3. Inspect the lockfile diff. With LAST and BLAST+ installed, run
   `uv run --no-sync python scripts/check.py --full --integration`. This checks
   version consistency, builds the selected version, and tests its distributions.
   Run both [worked examples](docs/README.md) when changing user-facing behavior.
4. Publish the checked commit to `master`. Administrators may push directly;
   other contributors must use a pull request and satisfy the branch's required
   checks and review conditions.
5. For an actual release, wait for CI to pass on that default-branch commit,
   then tag it as `vX.Y.Z` and describe the supported version and changes in the
   GitHub release.
   A version bump or passing check alone does not create a release or publish
   a package to PyPI. Resume development with a new `.devN` version and repeat
   the metadata/lockfile synchronization.

The maintained user documentation lives in [docs](docs/README.md). Update it,
the README, and the examples together when inputs, outputs, or CLI behavior
change. Keep any future Wiki navigation pointed at these files rather than
maintaining a second copy of the instructions.

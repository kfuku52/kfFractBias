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

Changes to synteny generation should also pass the real JCVI/LAST integration
test documented in the README. New parsing behavior must include malformed,
ambiguous, and duplicate-input tests. Never commit biological input data,
generated result directories, credentials, or local environments.

Synteny tests assert known pair sets, depth, diagonal boundaries, and retention
values, not just nonempty files. The fast suite includes a naive retention
oracle, README CLI parsing, subprocess hash-seed checks, multiprocess locking,
input mutation, source provenance, documentation examples/links, and
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
remains on Linux/Python 3.12. Distribution checks compare source fingerprints
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

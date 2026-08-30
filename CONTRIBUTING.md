# Contributing

Use Python 3.11 or newer and install the locked development environment:

```bash
uv sync --locked --extra test --extra plot
```

Before opening a pull request, run:

```bash
uv run --no-sync python scripts/check.py
```

For coverage, wheel/sdist build, and isolated distribution tests, use
`uv run --no-sync python scripts/check.py --full`. For changes to synteny,
install the locked `test` and `all` extras and LAST, then add `--integration`.
Use `--no-sync` while checking so an invocation does not remove optional
packages from an already prepared integration environment.

Changes to synteny generation should also pass the real JCVI/LAST integration
test documented in the README. New parsing behavior must include malformed,
ambiguous, and duplicate-input tests. Never commit biological input data,
generated result directories, credentials, or local environments.

Synteny tests assert known pair sets, depth, diagonal boundaries, and retention
values, not just nonempty files. The fast suite includes a naive retention
oracle, README CLI parsing, subprocess hash-seed checks, multiprocess locking,
input mutation, and commit/rollback failure injection. The optional JCVI core
tests do not require LAST. Use synthetic fixture generators for new tests.

Run `python scripts/benchmark.py --genes 10000 --queries 100` to measure dense
streaming tables. Compare equal Python/dependency environments and TSV hashes;
measure Python allocation peaks separately with `--trace-memory`, since the
tracer changes timing. Do not make noisy wall-clock thresholds a CI gate.

The quality job caches the same locked optional environment as integration.
Inspect setup/cache, aligner installation, and test-step durations separately
when assessing CI latency. A cache hit is not a guarantee of a fixed speedup.

Keep `pyproject.toml`, `CITATION.cff`, and `CHANGELOG.md` synchronized when
preparing a release. Releases are tagged as `vX.Y.Z` from the protected default
branch after all required checks pass.

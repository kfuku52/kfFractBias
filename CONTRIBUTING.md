# Contributing

Use Python 3.11 or newer and install the locked development environment:

```bash
uv sync --locked --extra test --extra plot
```

Before opening a pull request, run:

```bash
uv run --no-sync python -m pytest
uv run --no-sync ruff check src tests
uv run --no-sync ruff format --check src tests
uv run --no-sync mypy src
uv run --no-sync python -m build
```

Changes to synteny generation should also pass the real JCVI/LAST integration
test documented in the README. New parsing behavior must include malformed,
ambiguous, and duplicate-input tests. Never commit biological input data,
generated result directories, credentials, or local environments.

Keep `pyproject.toml`, `CITATION.cff`, and `CHANGELOG.md` synchronized when
preparing a release. Releases are tagged as `vX.Y.Z` from the protected default
branch after all required checks pass.

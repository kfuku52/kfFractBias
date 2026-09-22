---
name: check-cli-docs
description: Verify kfFractBias CLI and documentation changes against the parser and the existing temporary-directory example tests. Use for CLI flags, README/docs, or synthetic tutorial edits, not real-genome analysis or release validation.
---

# Check CLI documentation

Inputs: the changed CLI/document/example paths and the intended behavior.
Run from the repository root. Read the relevant diff and
[CONTRIBUTING](../../../CONTRIBUTING.md#choosing-checks) for setup and test tiers.

1. Identify affected commands and expected values in the maintained Markdown.
   Compare flags/defaults with `src/kffractbias/cli.py`; for changed numerical
   claims consult `docs/methods.md` and `docs/formats.md`.
2. Run the existing checks:

   ```bash
   uv run --no-sync python -m pytest -q -rs tests/test_documentation.py tests/test_cli.py
   ```

   `tests/documentation.py` discovers root, docs, and examples Markdown.
   It parses CLI examples without requiring placeholder input files. The two
   execution tests in `tests/test_documentation.py` copy/generate synthetic
   inputs under pytest temporary directories, run calculate/validate, and
   assert retention and mapping counts. Do not run the tutorial shell commands
   in the checkout as a substitute: their named output directories can persist
   or overwrite earlier results.
3. Review the coverage boundary. Shell fences must be bash/sh/shell/console.
   The parser checks recognized CLI invocations, not arbitrary shell programs;
   internal-link checks do not fetch external URLs. Hidden `.agents` Markdown
   is not discovered. Inspect its links manually when editing this skill.
   For compare/selfcompare, aligner, or plot changes, select the additional
   tests from CONTRIBUTING. Do not interpret parsed examples as executed
   alignments. Avoid weakening assertions to accommodate a regression.
4. Once edits are complete, run the standard check from CONTRIBUTING before
   push. No need to repeat a check already passed on the final diff.

Deliver a short verification record: changed documents/commands, exact check
command, pass/fail/skip results, and what was only inspected. Success requires
zero exit status and the relevant execution tests passing. If dependencies
are missing, use the documented locked environment; retain all extras in an
existing comparison environment. If required tools cannot be installed within
the task's scope, report that portion unverified and continue independent
checks. On failure, report the failing test and resolve the cause; never
suppress the failure or replace expected scientific results without evidence.

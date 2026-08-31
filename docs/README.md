# Documentation

These documents describe the maintained offline `kffractbias` CLI and are
versioned with its source code. Start with the repository [installation
instructions](../README.md). When reading an older result, use the documentation
at the commit that produced it; package versions alone did not distinguish all
historical untagged snapshots.

| Task | Reference |
| --- | --- |
| Check installation with BED and anchors | [Minimal example](../examples/minimal/README.md) |
| Run validate, compare, and selfcompare on synthetic CDS/GFF | [Annotation tutorial](../examples/annotations/README.md) |
| Understand ranks, denominators, windows, quota, and self retention | [Calculation rules](methods.md) |
| Prepare inputs or read TSV/JSON outputs | [Format reference](formats.md) |
| Diagnose installation, input, output, or alignment failures | [Troubleshooting](troubleshooting.md) |
| Update an old checkout or reproduce an old analysis | [Migration](migration.md) |
| Run development checks or prepare a release | [Contributing](../CONTRIBUTING.md) |
| Report a vulnerability privately | [Security policy](../SECURITY.md) |

The upstream CoGe Wiki documents the original web application and is historical
context for this fork. CoGe genome IDs, accounts, JWTs, and its web workflow are
not required here. The repository documents above are the maintained source of
instructions; a separate Wiki copy is not maintained.

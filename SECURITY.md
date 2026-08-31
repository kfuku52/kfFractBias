# Security policy

Before the first tagged release, security fixes are provided on the current
`master` development line. After a release is published, the latest tagged
release is supported. Historical untagged snapshots are not maintained separately.

Please [report suspected vulnerabilities privately on GitHub](https://github.com/kfuku52/kfFractBias/security/advisories/new)
rather than opening a public issue. Private vulnerability reporting is enabled
for this repository. Include the package version and commit, Python and tool
versions, a minimal reproduction, and the expected impact when possible. Keep
sensitive biological data and credentials out of public reports and logs.

The base `kffractbias` installation has no third-party runtime dependencies.
The optional `compare` extra installs JCVI and its scientific dependency tree;
use `uv sync --locked --extra all` for the reviewed dependency set and keep the
lock file updated. As reviewed on 2026-08-31, JCVI 1.6.6 pulls `diskcache 5.6.3`
through an unused `genomepy` dependency.
[CVE-2025-69872 / GHSA-w8v5-vhqr-4h9v](https://github.com/advisories/GHSA-w8v5-vhqr-4h9v)
has no fixed release and requires
an attacker to have write access to the cache directory; kfFractBias does not
import or create a DiskCache. CI allows only that advisory explicitly and
fails for any additional known vulnerability. Its pip-audit identifier is
`PYSEC-2026-2447`. Recheck the advisory when refreshing dependencies; the review
date is not a promise that upstream status will remain unchanged.

External FASTA, GFF, BED, and synteny files are untrusted inputs and should
not share a writable output path.

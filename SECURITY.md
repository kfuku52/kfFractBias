# Security policy

Security fixes are provided for the latest release. Please report suspected
vulnerabilities through GitHub's private vulnerability reporting for this
repository rather than a public issue. Include affected versions, a minimal
reproduction, and the expected impact when possible.

The base `kffractbias` installation has no third-party runtime dependencies.
The optional `compare` extra installs JCVI and its scientific dependency tree;
use `uv sync --locked --extra all` for the reviewed dependency set and keep the
lock file updated. JCVI 1.6.6 currently pulls `diskcache 5.6.3` through an
unused `genomepy` dependency. CVE-2025-69872 has no fixed release and requires
an attacker to have write access to the cache directory; kfFractBias does not
import or create a DiskCache. CI allows only that advisory explicitly and
fails for any additional known vulnerability.

External FASTA, GFF, BED, and synteny files are untrusted inputs and should
not share a writable output path.

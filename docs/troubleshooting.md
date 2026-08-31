# Troubleshooting

Start from the [installation instructions](../README.md) and keep the source
checkout separate from biological inputs and generated output directories.

| Symptom | Check or remedy |
| --- | --- |
| `kffractbias: command not found` after uv sync | Activate `.venv` with `source .venv/bin/activate`, or use `uv run --no-sync kffractbias ...` from the repository root. |
| No PyPI distribution for kffractbias | Install this source checkout. `python -m pip install '.[compare]'` must run at its root in the intended environment. |
| Missing matplotlib | Install the plot extra, or add `--no-plot` to analysis commands. |
| Missing JCVI or solver | Install compare/all extras in the environment used to execute the CLI. The solver backend is OR-Tools SCIP. |
| Missing `lastal` / `lastdb` | Install LAST and put both executables on PATH. They are not Python packages installed by uv/pip. |
| Missing `blastn` / `makeblastdb` | Install NCBI BLAST+ and put both executables on PATH, then use `--aligner blast`. |
| `--blast-task requires --aligner blast` | A task override only applies to BLAST. The default BLAST task is `blastn`. |
| Only part of CDS input maps to GFF | Check that FASTA IDs match the intended GFF feature/attribute and that the annotations belong to this FASTA. Feature/attribute overrides must be supplied together. |
| Multiple mapped isoforms of a locus | Use representative CDS input, `--isoform-policy longest`, or an intentional `all` policy. Do not relax the mapping fraction to address isoforms. |
| Unknown or overlapping pairwise gene IDs | Keep IDs consistent across FASTA/GFF/BED/anchors and disjoint between genomes. |
| No pairs remain after selection | Check sequence names and retained anchors. Selecting only unmatched sequences is an error even when unmatched output is requested. |
| windows.tsv contains only its header | Compare the window size with analyzed gene counts after filtering. Only complete windows are produced. |
| Missing chromosomes in the PNG | Open the PDF, which contains every page. PNG is only the first-page preview. |
| Huge TSV files despite low memory usage | Output remains dense over analyzed genes and included query sequences. See [row-count formulas](methods.md). |

For Ubuntu, CI uses `sudo apt-get install last-align ncbi-blast+` after updating
the package index. Other platforms may use a package manager or the upstream
tools' installers. Confirm availability in the same shell/environment as the
analysis with `command -v lastal`, `command -v lastdb`, `command -v blastn`, and
`command -v makeblastdb` as appropriate. Report `lastal --version` or
`blastn -version` when diagnosing alignment differences.

**Existing results and locks**

Successful analysis replaces result files with the same output directory and
prefix. A successful `--no-plot` rerun also removes that prefix's previous
figures. Choose a different prefix to preserve results. Compare/selfcompare
require `--force` to replace an existing synteny work directory, but force
never bypasses an active lock.

`Another analysis is writing the same output prefix` means another process
holds the lock. Wait for that process or use another prefix. Lock files remain
after release; their presence alone does not mean a process is running. Do not
delete or replace a lock file while a process may hold it, as that can split
locking across different files. Local POSIX filesystems are supported; Windows
and network filesystems are not supported or tested.

Inputs may not alias outputs or live in a replaceable synteny directory during
a comparison. Reuse retained BED/anchors with `calculate` and a new prefix,
rather than placing source FASTA/GFF inside work that `--force` will replace.
Allow disk space for original-input snapshots, newly generated results, and
old results kept until successful replacement.

**Failed comparisons**

Add `--keep-failed-work` to retain this run's unique `.PREFIX.staging-*`
directory. The CLI prints its location. Inspect `failure.json`, copied inputs,
and any partial `PREFIX.synteny/logs/` output. Successful comparisons retain
per-command logs and JSON records under the final synteny directory.

Ordinary failed replacements roll back to previous results. If rollback itself
fails, recovery files are retained automatically, even without the flag. Keep
those files until recovery is complete. A process killed abruptly or a power
failure does not have a multi-file atomicity guarantee; downstream readers
should wait for completion and verify output hashes in the final summary.

Diagnostic directories can contain genomic data. Remove only a known finished
run's diagnostic directory after it is no longer needed. For public bug
reports, prefer synthetic inputs and redact sensitive paths/content. Report
security problems through the [private reporting channel](../SECURITY.md).

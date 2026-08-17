# kfFractBias

`kfFractBias` is a Python 3 command-line program for calculating and plotting
gene-retention and fractionation-bias profiles. It is a maintained fork and
offline redesign of
[`SynMapFractBiasAnalysis`](https://github.com/bjoyce3/SynMapFractBiasAnalysis).

The original FractBias implementation was developed by Blake L. Joyce, Asher
Haug-Baltzell, Sean Davey, Matthew Bomhoff, James C. Schnable, and Eric Lyons.
The original Python 2 and notebook implementations remain under
`Code_FractBias/` for provenance.

## What is different

- Python 3 package and the `kffractbias` executable
- no CoGe account, genome ID, API, or JWT requirement
- calculation from either JCVI anchors or legacy SynMap/DAGCHAINER output
- end-to-end local comparison from target/query CDS FASTA and GFF annotations
- local JCVI MCscan and QUOTA-ALIGN execution
- explicit target:query syntenic-depth quota
- long-form gene and sliding-window TSV outputs
- JSON provenance with input hashes and parameters
- non-interactive PDF and PNG plots

After the software and its dependencies are installed, analysis does not need
network access.

## Installation

```bash
python -m pip install .
kffractbias version
```

The `compare` command requires JCVI and one of its supported local aligners.
The default `last` aligner requires `lastal` and `lastdb`. JCVI QUOTA-ALIGN
also requires its mixed-integer solver dependency.

## Commands

### Calculate from precomputed synteny

`calculate` accepts a JCVI `.anchors` file or SynMap/DAGCHAINER output. The
target and query BED identifiers must match the synteny identifiers.

```bash
kffractbias calculate \
  --synteny target.query.lifted.1x2.anchors \
  --format jcvi \
  --target-bed target.bed \
  --query-bed query.bed \
  --target-name Sorghum_bicolor \
  --query-name Zea_mays \
  --window-size 100 \
  --denominator all \
  --output-dir results \
  --prefix sorghum_maize
```

### Compare two annotations end-to-end

`compare` detects matching GFF feature/attribute identifiers, prepares JCVI
BED and CDS inputs, runs MCscan and QUOTA-ALIGN, and then calculates the
fractionation profile.

```bash
kffractbias compare \
  --target-cds target.cds.fa \
  --target-gff target.gff3 \
  --query-cds query.cds.fa \
  --query-gff query.gff3 \
  --target-name Sorghum_bicolor \
  --query-name Zea_mays \
  --quota 1:2 \
  --window-size 100 \
  --cpus 8 \
  --output-dir results \
  --prefix sorghum_maize
```

The quota is always explicit: `1:2` means one expected target region for two
expected query regions. Reversing target and query also requires reversing the
quota.

### Validate annotations

```bash
kffractbias validate \
  --target-cds target.cds.fa --target-gff target.gff3 \
  --query-cds query.cds.fa --query-gff query.gff3
```

Use `--target-feature` with `--target-attribute`, and the corresponding query
options, when automatic GFF identifier detection is not appropriate.

Other commands are `kffractbias formats` and `kffractbias version`.

## Outputs

For a prefix such as `sorghum_maize`, kfFractBias writes:

- `sorghum_maize.genes.tsv`: per-target-gene retention and matching query IDs
- `sorghum_maize.windows.tsv`: sliding-window retention fractions and percent
- `sorghum_maize.summary.json`: parameters, counts, input hashes, and outputs
- `sorghum_maize.plot.pdf` and `.plot.png`: chromosome-wise profiles
- `sorghum_maize.synteny/`: retained JCVI working data from `compare`

`--denominator all` uses every target gene in each window. `--denominator
syntenic` first removes target genes without any retained query match, matching
the two denominator choices exposed by the original FractBias implementation.

## Development

```bash
python -m pip install -e '.[test]'
pytest
```

## Citation

If you use kfFractBias, cite the original FractBias publication:

> Joyce BL, Haug-Baltzell A, Davey S, Bomhoff M, Schnable JC, Lyons E.
> FractBias: a graphical tool for assessing fractionation bias following
> polyploidy. Bioinformatics. 2017;33(4):552–554.
> <https://doi.org/10.1093/bioinformatics/btw666>

## License

MIT. The original copyright and license are retained in `MIT License`.

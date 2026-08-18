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
- within-genome self-synteny retention analysis with native JCVI self filtering
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

### Compare a genome to itself

`selfcompare` measures within-genome retention across non-redundant self-synteny
blocks. It runs JCVI with the same BED/CDS prefix on both axes so JCVI removes
identity hits, mirrored pairs, tandem-proximal hits, and the near-diagonal
intrachromosomal region before symmetric QUOTA-ALIGN screening.

```bash
kffractbias selfcompare \
  --cds polyploid.cds.fa \
  --gff polyploid.gff3 \
  --name Polyploid_species \
  --depth 1 \
  --window-size 100 \
  --cpus 8 \
  --output-dir results \
  --prefix polyploid_self
```

`--depth N` applies a symmetric `N:N` maximum block depth. It must be chosen
from the expected homeologous structure rather than inferred automatically.
The default `--self-hit-percent 98` is JCVI's near-self identity cutoff and can
be changed for unusually recent polyploidy. JCVI also excludes
intrachromosomal anchors within 300 gene ranks of the self diagonal by default;
use `--diagonal-bound` to change that threshold for genomes with short
chromosomes or for a deliberately different diagonal definition.

The result is explicitly recorded as `self_synteny_retention`. It measures
retention asymmetry conditional on genes present in the annotated genome; it
cannot observe ancestral loci for which every descendant copy was lost and is
therefore not equivalent to outgroup-based fractionation bias from `compare`.

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

- `sorghum_maize.genes.tsv`: retention and matching query IDs for each target
  gene and query sequence combination
- `sorghum_maize.windows.tsv`: sliding-window retention fractions and percent
- `sorghum_maize.summary.json`: parameters, counts, input hashes, and outputs
- `sorghum_maize.plot.pdf` and `.plot.png`: chromosome-wise profiles
- `sorghum_maize.synteny/`: retained JCVI working data from `compare`

`selfcompare` writes the same output set and records identity/mirror filtering,
interchromosomal and intrachromosomal pair counts, symmetric depth, and the
interpretation limitation in the JSON summary.

`--denominator all` uses every target gene in each window. `--denominator
syntenic` first removes target genes without any retained query match, matching
the two denominator choices exposed by the original FractBias implementation.

## Development

```bash
python -m pip install -e '.[test]'
python -m pytest
python -m ruff check src tests
```

## Citation

If you use kfFractBias, cite the original FractBias publication:

> Joyce BL, Haug-Baltzell A, Davey S, Bomhoff M, Schnable JC, Lyons E.
> FractBias: a graphical tool for assessing fractionation bias following
> polyploidy. Bioinformatics. 2017;33(4):552–554.
> <https://doi.org/10.1093/bioinformatics/btw666>

## License

MIT. The original copyright and license are retained in `MIT License`.

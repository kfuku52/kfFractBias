from pathlib import Path

import pytest

from kffractbias.io import Gene, parse_synteny_pairs, read_bed, write_bed
from kffractbias.selfscan import scan_self, select_quota


def alignment(path: Path, pairs: list[tuple[str, str]]) -> Path:
    path.write_text(
        "".join(f"{a}\t{b}\t95\t300\t15\t0\t1\t300\t1\t300\t1e-40\t200\n" for a, b in pairs)
    )
    return path


def test_solver_enforces_depth_across_both_axes():
    pytest.importorskip("jcvi")
    genes = tuple(Gene(chrom, i, i + 1, f"{chrom}{i}") for chrom in "ABC" for i in range(8))
    blocks = [[(i + a, i + b, 200) for i in range(8)] for a, b in ((0, 8), (0, 16), (8, 16))]
    assert len(select_quota(blocks, genes, 1)) == 1
    assert select_quota(blocks, genes, 2) == [0, 1, 2]


@pytest.mark.parametrize("distance,expected", [(299, 0), (300, 4)])
def test_self_diagonal_boundary_and_liftover(tmp_path, distance, expected):
    pytest.importorskip("jcvi")
    genes = tuple(Gene("chr1", i, i + 1, f"g{i}") for i in range(distance + 4))
    bed = tmp_path / "self.bed"
    write_bed(genes, bed)
    pairs = [(f"g{i}", f"g{i + distance}") for i in range(4)]
    # Include identity, mirrors, and a near-diagonal liftover candidate.
    blast = alignment(tmp_path / "hits.last", pairs + [(b, a) for a, b in pairs] + [("g0", "g0")])
    lift = alignment(tmp_path / "lift.last", pairs + [("g3", "g302")])
    output = tmp_path / "self.self.lifted.1x1.anchors"
    if expected:
        scan_self(blast, lift, bed, output, bound=300, depth=1)
        assert set(
            parse_synteny_pairs(
                output,
                "jcvi",
                {g.gene_id for g in genes},
                {g.gene_id for g in genes},
                allow_ambiguous_orientation=True,
            ).pairs
        ) == set(pairs)
    else:
        with pytest.raises(ValueError, match="No self-synteny blocks"):
            scan_self(blast, lift, bed, output, bound=300, depth=1)


@pytest.mark.parametrize("padding_chrom", [None, "chr0"])
def test_interchromosomal_pairs_ignore_unrelated_chromosomes(tmp_path, padding_chrom):
    pytest.importorskip("jcvi")
    genes = [Gene(chrom, i, i + 1, f"{chrom}{i}") for chrom in ("chrA", "chrB") for i in range(8)]
    if padding_chrom:
        genes += [Gene(padding_chrom, i, i + 1, f"padding{i}") for i in range(320)]
    bed = tmp_path / "self.bed"
    write_bed(reversed(genes), bed)
    expected = [(f"chrA{i}", f"chrB{i}") for i in range(8)]
    blast = alignment(tmp_path / "hits.last", expected + [(b, a) for a, b in expected])
    output = tmp_path / "self.self.lifted.1x1.anchors"
    scan_self(blast, blast, bed, output, bound=300, depth=1)
    identifiers = {gene.gene_id for gene in read_bed(bed)}
    assert set(
        parse_synteny_pairs(
            output, "jcvi", identifiers, identifiers, allow_ambiguous_orientation=True
        ).pairs
    ) == set(expected)

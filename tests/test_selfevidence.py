import csv
import json
from collections import Counter

import pytest

from kffractbias.io import Gene, sha256_file
from kffractbias.selfevidence import write_self_evidence
from kffractbias.selfscan import write_blocks


def evidence(tmp_path, blocks, genes):
    write_blocks(tmp_path / "self.self.anchors", blocks, genes)
    write_blocks(tmp_path / "self.self.lifted.anchors", blocks, genes)
    paths = write_self_evidence(blocks, genes, tmp_path)
    with paths["depth"].open() as handle:
        depth = list(csv.DictReader(handle, delimiter="\t"))
    with paths["blocks"].open() as handle:
        block_rows = list(csv.DictReader(handle, delimiter="\t"))
    return paths, depth, block_rows, json.loads(paths["summary"].read_text())


def test_raw_depth_coverage_zero_rows_and_hashes(tmp_path):
    genes = tuple(
        Gene(chrom, i * 10, i * 10 + 5, f"{chrom}{i}") for chrom in "ABC" for i in range(5)
    )
    blocks = [[(0, 5, 100), (2, 7, 90)], [(1, 11, 100), (3, 13, 80)]]
    paths, depth, rows, summary = evidence(tmp_path, blocks, genes)
    assert [int(row["block_arm_depth"]) for row in depth[:5]] == [1, 2, 2, 1, 0]
    assert [int(row["num_anchor_partners"]) for row in depth[:5]] == [1, 1, 1, 1, 0]
    assert int(depth[2]["num_partner_chromosomes"]) == 2
    assert rows[0]["b_start_rank"] == "1"
    assert rows[0]["b_end_rank"] == "3"
    assert rows[0]["b_span_genes"] == "3"
    assert summary["num_genes"] == 15
    assert summary["num_span_covered_genes"] == 10
    assert summary["num_anchor_genes"] == 8
    assert summary["gene_span_coverage"] == pytest.approx(10 / 15)
    assert summary["block_arm_depth_distribution"] == {"0": 5, "1": 8, "2": 2}
    assert "not ploidy" in summary["interpretation"]
    for key, digest in summary["output_sha256"].items():
        assert digest == sha256_file(paths[key])


def test_overlapping_arms_count_twice_and_duplicate_anchors_once(tmp_path):
    genes = tuple(Gene("chr1", i, i + 1, f"g{i}") for i in range(8))
    blocks = [[(0, 2, 100), (4, 6, 100), (4, 6, 100)]]
    _, depth, rows, _ = evidence(tmp_path, blocks, genes)
    assert [int(row["block_arm_depth"]) for row in depth] == [1, 1, 2, 2, 2, 1, 1, 0]
    assert rows[0]["num_anchor_rows"] == "3"
    assert rows[0]["num_unique_anchor_pairs"] == "2"
    assert depth[4]["num_anchor_partners"] == "1"


def test_depth_matches_direct_coverage_with_overlapping_partner_intervals(tmp_path):
    genes = tuple(Gene(chrom, i, i + 1, f"{chrom}{i}") for chrom in "ABC" for i in range(6))
    blocks = [
        [(0, 6, 100), (2, 8, 90)],
        [(1, 7, 100), (3, 9, 90)],
        [(3, 12, 100), (5, 14, 90)],
        [(4, 10, 80), (4, 10, 80)],
    ]
    _, rows, _, summary = evidence(tmp_path, blocks, genes)
    depths = [0] * len(genes)
    chromosomes = [set() for _ in genes]
    partners = [set() for _ in genes]
    for block in blocks:
        for left, right, _ in block:
            partners[left].add(right)
            partners[right].add(left)
        for axis in (0, 1):
            ranks = [point[axis] for point in block]
            opposite_chromosome = genes[block[0][1 - axis]].seqid
            for rank in range(min(ranks), max(ranks) + 1):
                depths[rank] += 1
                chromosomes[rank].add(opposite_chromosome)
    assert [int(row["block_arm_depth"]) for row in rows] == depths
    assert [int(row["num_partner_chromosomes"]) for row in rows] == list(map(len, chromosomes))
    assert [int(row["num_anchor_partners"]) for row in rows] == list(map(len, partners))
    assert summary["num_span_covered_genes"] == sum(depth > 0 for depth in depths)
    assert summary["block_arm_depth_distribution"] == {
        str(depth): count for depth, count in Counter(depths).items()
    }


def test_empty_blocks_preserve_all_zero_rows(tmp_path):
    gene_count = 3
    genes = tuple(Gene("A", i, i + 1, f"g{i}") for i in range(gene_count))
    _, rows, block_rows, summary = evidence(tmp_path, [], genes)
    assert len(rows) == gene_count
    assert block_rows == []
    assert all(row["block_arm_depth"] == row["num_partner_chromosomes"] == "0" for row in rows)
    assert summary["gene_span_coverage"] == 0
    assert summary["block_arm_depth_distribution"] == {"0": gene_count}

import csv
import json

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


def test_invalid_block_coordinates_fail(tmp_path):
    genes = (Gene("A", 0, 1, "a"), Gene("B", 0, 1, "b"))
    with pytest.raises(ValueError, match="outside"):
        write_self_evidence([[(0, 2, 1)]], genes, tmp_path)
    with pytest.raises(ValueError, match="within chromosomes"):
        write_self_evidence([[(0, 1, 1), (1, 0, 1)]], genes, tmp_path)

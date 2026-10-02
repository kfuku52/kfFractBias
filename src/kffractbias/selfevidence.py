"""Unquota self-block coverage audits; these are not ploidy estimates."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from .io import Gene, sha256_file

Point = tuple[int, int, int]


def write_self_evidence(
    blocks: list[list[Point]], genes: tuple[Gene, ...], directory: Path
) -> dict[str, Path]:
    paths = {
        "scan_anchors": directory / "self.self.anchors",
        "lifted_anchors": directory / "self.self.lifted.anchors",
        "blocks": directory / "self.self.raw.blocks.tsv",
        "depth": directory / "self.self.raw.depth.tsv",
        "summary": directory / "self.self.raw.summary.json",
    }
    local_ranks: dict[int, int] = {}
    chromosome_counts: Counter[str] = Counter()
    for index, gene in enumerate(genes):
        chromosome_counts[gene.seqid] += 1
        local_ranks[index] = chromosome_counts[gene.seqid]
    spans: Counter[int] = Counter()
    partner_chromosomes: dict[int, set[str]] = defaultdict(set)
    anchor_partners: dict[int, set[int]] = defaultdict(set)
    block_rows = []
    for index, block in enumerate(blocks, 1):
        if not block:
            raise ValueError("Self-evidence blocks must be nonempty.")
        axes = [[point[axis] for point in block] for axis in (0, 1)]
        if any(rank < 0 or rank >= len(genes) for axis in axes for rank in axis):
            raise ValueError("Self-evidence anchor rank is outside the prepared BED.")
        if any(len({genes[rank].seqid for rank in axis}) != 1 for axis in axes):
            raise ValueError("Self-evidence block arms must remain within chromosomes.")
        for left, right, _score in block:
            anchor_partners[left].add(right)
            anchor_partners[right].add(left)
        row: dict[str, object] = {
            "block_id": f"block{index:06d}",
            "num_anchor_rows": len(block),
            "num_unique_anchor_pairs": len({(left, right) for left, right, _ in block}),
        }
        for arm, ranks, opposite in zip(("a", "b"), axes, reversed(axes), strict=True):
            low, high = min(ranks), max(ranks)
            chromosome = genes[low].seqid
            partner = genes[opposite[0]].seqid
            for rank in range(low, high + 1):
                if genes[rank].seqid != chromosome:
                    raise ValueError("Prepared BED chromosomes must be contiguous in rank order.")
                spans[rank] += 1
                partner_chromosomes[rank].add(partner)
            row.update(
                {
                    f"{arm}_seqid": chromosome,
                    f"{arm}_start_rank": local_ranks[low],
                    f"{arm}_end_rank": local_ranks[high],
                    f"{arm}_start_gene": genes[low].gene_id,
                    f"{arm}_end_gene": genes[high].gene_id,
                    f"{arm}_start": min(genes[rank].start for rank in ranks),
                    f"{arm}_end": max(genes[rank].end for rank in ranks),
                    f"{arm}_num_anchor_genes": len(set(ranks)),
                    f"{arm}_span_genes": high - low + 1,
                }
            )
        block_rows.append(row)
    block_columns = ["block_id", "num_anchor_rows", "num_unique_anchor_pairs"] + [
        f"{arm}_{field}"
        for arm in ("a", "b")
        for field in (
            "seqid",
            "start_rank",
            "end_rank",
            "start_gene",
            "end_gene",
            "start",
            "end",
            "num_anchor_genes",
            "span_genes",
        )
    ]
    with paths["blocks"].open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=block_columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(block_rows)
    depth_columns = [
        "seqid",
        "gene_id",
        "gene_rank",
        "start",
        "end",
        "block_arm_depth",
        "num_anchor_partners",
        "num_partner_chromosomes",
    ]
    with paths["depth"].open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=depth_columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for rank, gene in enumerate(genes):
            writer.writerow(
                {
                    "seqid": gene.seqid,
                    "gene_id": gene.gene_id,
                    "gene_rank": local_ranks[rank],
                    "start": gene.start,
                    "end": gene.end,
                    "block_arm_depth": spans[rank],
                    "num_anchor_partners": len(anchor_partners.get(rank, set())),
                    "num_partner_chromosomes": len(partner_chromosomes.get(rank, set())),
                }
            )
    distribution = Counter(spans[rank] for rank in range(len(genes)))
    summary = {
        "schema_version": 1,
        "stage": "prequota-lifted-self-blocks",
        "num_genes": len(genes),
        "num_blocks": len(blocks),
        "num_span_covered_genes": sum(spans[rank] > 0 for rank in range(len(genes))),
        "num_anchor_genes": len(anchor_partners),
        "gene_span_coverage": sum(spans[rank] > 0 for rank in range(len(genes))) / len(genes)
        if genes
        else 0,
        "anchor_gene_coverage": len(anchor_partners) / len(genes) if genes else 0,
        "block_arm_depth_distribution": {
            str(depth): count for depth, count in sorted(distribution.items())
        },
        "outputs": {key: path.name for key, path in paths.items() if key != "summary"},
        "output_sha256": {
            key: sha256_file(path) for key, path in paths.items() if key != "summary"
        },
        "interpretation": "Overlapping unquota block-arm spans on the full prepared BED, not ploidy, homoeology, ancestral copy number or a WGD test. Redundant/local blocks can increase depth. Both arms count separately where they overlap.",
    }
    paths["summary"].write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )
    return paths

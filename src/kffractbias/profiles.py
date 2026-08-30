"""Typed, streaming retention tables, independent of files and run state."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from typing import TypedDict

from .io import Gene, natural_key


class GeneRow(TypedDict):
    target_seqid: str
    target_gene: str
    target_rank: int
    query_seqid: str
    retained: int
    query_genes: str


class WindowRow(TypedDict):
    target_seqid: str
    query_seqid: str
    window_index: int
    start_rank: int
    end_rank: int
    start_gene: str
    end_gene: str
    retained_count: int
    window_size: int
    retention_fraction: str
    retention_percent: str


GENE_FIELDS = tuple(GeneRow.__annotations__)
WINDOW_FIELDS = tuple(WindowRow.__annotations__)


@dataclass(frozen=True)
class RetentionProfile:
    target_by_seqid: dict[str, list[Gene]]
    query_seqids: list[str]
    mappings: dict[str, dict[str, set[str]]]
    window_size: int
    step_size: int

    def gene_rows(self) -> Iterator[GeneRow]:
        for seqid in sorted(self.target_by_seqid, key=natural_key):
            for rank, gene in enumerate(self.target_by_seqid[seqid], start=1):
                matches = self.mappings.get(gene.gene_id, {})
                for query_seqid in self.query_seqids:
                    query_ids = matches.get(query_seqid, ())
                    yield {
                        "target_seqid": seqid,
                        "target_gene": gene.gene_id,
                        "target_rank": rank,
                        "query_seqid": query_seqid,
                        "retained": int(bool(query_ids)),
                        "query_genes": ";".join(sorted(query_ids, key=natural_key)),
                    }

    def window_rows(self) -> Iterator[WindowRow]:
        for seqid in sorted(self.target_by_seqid, key=natural_key):
            genes = self.target_by_seqid[seqid]
            if len(genes) < self.window_size:
                continue
            for query_seqid in self.query_seqids:
                cumulative = [0]
                for gene in genes:
                    cumulative.append(
                        cumulative[-1]
                        + int(bool(self.mappings.get(gene.gene_id, {}).get(query_seqid)))
                    )
                for start in range(0, len(genes) - self.window_size + 1, self.step_size):
                    end = start + self.window_size
                    retained = cumulative[end] - cumulative[start]
                    fraction = retained / self.window_size
                    yield {
                        "target_seqid": seqid,
                        "query_seqid": query_seqid,
                        "window_index": start // self.step_size + 1,
                        "start_rank": start + 1,
                        "end_rank": end,
                        "start_gene": genes[start].gene_id,
                        "end_gene": genes[end - 1].gene_id,
                        "retained_count": retained,
                        "window_size": self.window_size,
                        "retention_fraction": f"{fraction:.10g}",
                        "retention_percent": f"{fraction * 100:.10g}",
                    }

    def row_counts(self) -> tuple[int, int]:
        queries = len(self.query_seqids)
        genes = sum(map(len, self.target_by_seqid.values())) * queries
        windows = (
            sum(
                max(0, (len(chrom) - self.window_size) // self.step_size + 1)
                for chrom in self.target_by_seqid.values()
            )
            * queries
        )
        return genes, windows

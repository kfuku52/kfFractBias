"""Chromosome-aware adapter for JCVI's self scan and QUOTA-ALIGN solver.

JCVI 1.6 applies the self diagonal bound to concatenated ranks even across
chromosomes. Keep its chaining and liftover algorithms, but filter hits before
either operation. Quota constraints share one genome axis, so a region cannot
evade its depth limit by appearing on different sides of different blocks.
Optional JCVI imports stay inside this subprocess entry point and its helpers.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from .io import Gene, natural_key, read_bed

Point = tuple[int, int, int]
Block = list[Point]


def align_self(aligner: str, cpus: int, cscore: float, self_hit_percent: float) -> None:
    """Use JCVI's alignment/filter stages without running its incorrect scan."""
    from jcvi.apps.align import last
    from jcvi.compara.blastfilter import main as filter_synteny_hits
    from jcvi.formats.blast import filter as filter_identity_hits
    from jcvi.formats.blast import filtered_blastfile_name

    if (
        aligner not in {"last", "blast"}
        or cpus < 1
        or not 0 < cscore <= 1
        or not 0 < self_hit_percent <= 100
    ):
        raise ValueError("Invalid self-alignment parameters")
    raw = "self.self.last"
    if aligner == "last":
        last(["self.cds", "self.cds", f"--cpus={cpus}"], "nucl")
    elif not Path(raw).is_file():
        raise ValueError("BLAST output must be prepared before self filtering")
    inverse = filtered_blastfile_name(raw, self_hit_percent, 0, inverse=True)
    filter_identity_hits(
        [raw, "--hitlen=0", f"--pctid={self_hit_percent}", "--inverse", "--noself"]
    )
    filter_synteny_hits([inverse, f"--cscore={cscore}", "--tandem_Nmax=10", "--no_strip_names"])


def read_self_hits(path: Path, genes: tuple[Gene, ...], bound: int) -> dict[tuple[str, str], Block]:
    from jcvi.formats.blast import Blast

    order = {gene.gene_id: index for index, gene in enumerate(genes)}
    scores: dict[tuple[int, int], int] = {}
    for hit in Blast(str(path)):
        if hit.query not in order or hit.subject not in order:
            continue
        left, right = sorted((order[hit.query], order[hit.subject]))
        if left == right:
            continue
        if genes[left].seqid == genes[right].seqid and right - left < bound:
            continue
        pair = (left, right)
        scores[pair] = max(scores.get(pair, 0), int(hit.score))
    grouped: dict[tuple[str, str], Block] = defaultdict(list)
    for (left, right), score in sorted(scores.items()):
        grouped[(genes[left].seqid, genes[right].seqid)].append((left, right, score))
    return grouped


def quota_constraints(
    blocks: list[Block], genes: tuple[Gene, ...], depth: int
) -> list[tuple[int, ...]]:
    """Sweep inclusive block intervals on both arms of the same genome."""
    events: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for block_id, block in enumerate(blocks):
        for axis in (0, 1):
            ranks = [point[axis] for point in block]
            seqid = genes[ranks[0]].seqid
            events[seqid].extend(((min(ranks), 1, block_id), (max(ranks) + 1, -1, block_id)))
    constraints: set[tuple[int, ...]] = set()
    for seqid in sorted(events, key=natural_key):
        active: dict[int, int] = defaultdict(int)
        for _rank, change, block_id in sorted(events[seqid]):
            active[block_id] += change
            if not active[block_id]:
                del active[block_id]
            if len(active) > depth:
                constraints.add(tuple(sorted(active)))
    return sorted(constraints)


def select_quota(blocks: list[Block], genes: tuple[Gene, ...], depth: int) -> list[int]:
    from jcvi.algorithms.lpsolve import MIPDataModel
    from ortools.linear_solver import pywraplp

    if not blocks:
        return []
    constraints = quota_constraints(blocks, genes, depth)
    scores = [min(len({p[0] for p in block}), len({p[1] for p in block})) for block in blocks]
    model = MIPDataModel(
        [dict.fromkeys(constraint, 1) for constraint in constraints],
        [depth] * len(constraints),
        scores,
        len(blocks),
        len(constraints),
    )
    solver, variables = model.create_solver()
    if solver.Solve() != pywraplp.Solver.OPTIMAL:
        raise RuntimeError("Self QUOTA-ALIGN did not find an optimal solution")
    return [index for index in range(len(blocks)) if variables[index].solution_value() > 0.5]


def write_blocks(path: Path, blocks: list[Block], genes: tuple[Gene, ...]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for block in blocks:
            handle.write("###\n")
            for left, right, score in sorted(block):
                handle.write(f"{genes[left].gene_id}\t{genes[right].gene_id}\t{score}\n")


def scan_self(
    filtered: Path, liftover: Path, bed: Path, output: Path, *, bound: int, depth: int
) -> None:
    from jcvi.compara.synteny import synteny_liftover, synteny_scan

    if bound < 1 or depth < 1:
        raise ValueError("Self diagonal bound and depth must be positive")
    genes = read_bed(bed)
    hits = read_self_hits(filtered, genes, bound)
    blocks: list[Block] = []
    for chroms in sorted(hits, key=lambda pair: (natural_key(pair[0]), natural_key(pair[1]))):
        # Diagonal, identity and mirror filtering has already been done on hits.
        blocks.extend(synteny_scan(hits[chroms], 20, 20, 4, is_self=False))
    if not blocks:
        raise ValueError("No self-synteny blocks survived the requested diagonal bound")
    write_blocks(output.parent / "self.self.anchors", blocks, genes)

    all_hits = read_self_hits(liftover, genes, bound)
    anchors: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    membership: dict[tuple[int, int], int] = {}
    for block_id, block in enumerate(blocks):
        for left, right, _score in block:
            anchors[(genes[left].seqid, genes[right].seqid)].append((left, right))
            membership[(left, right)] = block_id
    for chroms, points in anchors.items():
        if not all_hits.get(chroms):
            continue
        for point, nearest in synteny_liftover(all_hits[chroms], points, 10):
            left, right, score = map(int, point)
            blocks[membership[nearest]].append((left, right, score))

    write_blocks(output.parent / "self.self.lifted.anchors", blocks, genes)
    selected = select_quota(blocks, genes, depth)
    if not selected:
        raise ValueError("No self-synteny blocks survived QUOTA-ALIGN")
    write_blocks(output, [blocks[index] for index in selected], genes)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    align = commands.add_parser("align")
    align.add_argument("--aligner", choices=("last", "blast"), required=True)
    align.add_argument("--cpus", type=int, required=True)
    align.add_argument("--cscore", type=float, required=True)
    align.add_argument("--self-hit-percent", type=float, required=True)
    scan = commands.add_parser("scan")
    scan.add_argument("filtered", type=Path)
    scan.add_argument("liftover", type=Path)
    scan.add_argument("bed", type=Path)
    scan.add_argument("output", type=Path)
    scan.add_argument("--diagonal-bound", type=int, required=True)
    scan.add_argument("--depth", type=int, required=True)
    args = parser.parse_args()
    if args.command == "align":
        align_self(args.aligner, args.cpus, args.cscore, args.self_hit_percent)
        return
    scan_self(
        args.filtered,
        args.liftover,
        args.bed,
        args.output,
        bound=args.diagonal_bound,
        depth=args.depth,
    )


if __name__ == "__main__":
    main()

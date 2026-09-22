import gzip
from pathlib import Path

import pytest

from kffractbias.io import (
    annotation_to_genes,
    read_bed,
    read_fasta_ids,
    read_jcvi_pairs,
    read_synmap_pairs,
)
from kffractbias.jcvi import prepare_genome


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_annotation_mapping_detects_feature_and_attribute(tmp_path):
    fasta = write(tmp_path / "genes.fa", ">tx1\nATG\n>tx2\nATG\n")
    gff = write(
        tmp_path / "genes.gff3",
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1\t12\t.\t+\t.\tID=tx%31;Parent=g1\n"
        "chr1\ttest\tmRNA\t20\t30\t.\t-\t.\tID=tx2;Parent=g1\n"
        "##FASTA\n>chr1\nATG\n",
    )
    mapping = annotation_to_genes(gff, read_fasta_ids(fasta))
    assert mapping.feature == "mRNA"
    assert mapping.attribute == "ID"
    assert mapping.matched_gene_count == 2
    assert [(gene.gene_id, gene.start, gene.end) for gene in mapping.genes] == [
        ("tx1", 0, 12),
        ("tx2", 19, 30),
    ]


def test_duplicate_fasta_identifiers_are_rejected(tmp_path):
    fasta = write(tmp_path / "genes.fa", ">tx1\nATG\n>tx1 duplicate\nATG\n")
    with pytest.raises(ValueError, match="Duplicate FASTA identifier"):
        read_fasta_ids(fasta)


def test_synmap_parser_accepts_coge_numeric_identifiers(tmp_path):
    synmap = write(
        tmp_path / "synmap.txt",
        "a1_chr1\tchr1||1||2||left_name||1||CDS||101||1||99\t1\t2\t"
        "b2_chrA\tchrA||3||4||right_name||1||CDS||202||1||99\t3\t4\t1e-20\t10\n",
    )
    assert read_synmap_pairs(synmap, {"202"}, {"101"}) == (("202", "101"),)


@pytest.mark.parametrize(
    ("text", "message"),
    (
        ("only-one-column\n", "Expected at least 2 JCVI columns"),
        ("unknown\tq1\t10\n", "does not match one target and one query"),
    ),
)
def test_jcvi_parser_rejects_every_invalid_data_row(tmp_path, text, message):
    anchors = write(tmp_path / "pairs.anchors", f"t1\tq1\t10\n{text}")
    with pytest.raises(ValueError, match=message):
        read_jcvi_pairs(anchors, {"t1"}, {"q1"})


def test_jcvi_parser_rejects_ambiguous_orientation(tmp_path):
    anchors = write(tmp_path / "pairs.anchors", "g1\tg2\t10\n")
    with pytest.raises(ValueError, match="Ambiguous JCVI pair orientation"):
        read_jcvi_pairs(anchors, {"g1", "g2"}, {"g1", "g2"})


@pytest.mark.parametrize(
    "text,line,message",
    [
        (">\nATG\n", 1, "Empty FASTA identifier"),
        ("ATG\n>x\nATG\n", 1, "before first FASTA header"),
        (">x\n>y\nATG\n", 1, "Empty FASTA sequence"),
        (">x\n", 1, "Empty FASTA sequence"),
        (">x\nATX\n", 2, "Invalid nucleotide"),
    ],
)
def test_malformed_fasta_has_path_and_line(tmp_path, text, line, message):
    fasta = write(tmp_path / "invalid.fa", text)
    with pytest.raises(ValueError, match=message) as caught:
        read_fasta_ids(fasta)
    assert f"{fasta}:{line}" in str(caught.value)


def test_wrapped_compressed_iupac_fasta(tmp_path):
    fasta = tmp_path / "cds.fa.gz"
    with gzip.open(fasta, "wt") as handle:
        handle.write(">a description\nacgt\n nryswkmbdhvu\n\n>b\nATG\n")
    assert read_fasta_ids(fasta) == {"a", "b"}


@pytest.mark.parametrize("format", ["gff", "gtf", "cds"])
def test_isoform_policies_and_mapping_counts(tmp_path, format):
    fasta = write(tmp_path / "cds.fa", ">t1\nATG\n>t2\nATGATG\n>t3\nATGATG\n>t4\nATGATG\n")
    lines = []
    for i in range(1, 5):
        locus = "g1" if i <= 2 else "g2"
        if format == "gtf":
            lines.append(
                f'chr1\ttest\ttranscript\t{i * 10}\t{i * 10 + 5}\t.\t+\t.\tgene_id "{locus}"; transcript_id "t{i}";\n'
            )
        else:
            lines.append(
                f"chr1\ttest\tmRNA\t{i * 10}\t{i * 10 + 5}\t.\t+\t.\tID=t{i};Parent={locus}\n"
            )
            if format == "cds":
                lines.append(f"chr1\ttest\tCDS\t{i * 10}\t{i * 10 + 5}\t.\t+\t0\tParent=t{i}\n")
    gff = write(tmp_path / "genes.gff", "".join(lines))
    overrides = {"feature": "CDS", "attribute": "Parent"} if format == "cds" else {}
    with pytest.raises(ValueError, match="Multiple CDS isoforms"):
        prepare_genome("x", fasta, gff, tmp_path, **overrides)
    longest = prepare_genome("x", fasta, gff, tmp_path, isoform_policy="longest", **overrides)
    assert [g.gene_id for g in longest.mapping.genes] == ["t2", "t3"]
    assert read_fasta_ids(longest.cds_path) == {"t2", "t3"}
    assert longest.mapping.matched_gene_count == 4
    assert longest.mapping.collapsed_isoform_count == 2
    assert longest.mapping.metadata()["selected_gene_count"] == 2
    all_isoforms = prepare_genome("y", fasta, gff, tmp_path, isoform_policy="all", **overrides)
    assert len(all_isoforms.mapping.genes) == 4
    assert all_isoforms.mapping.metadata()["counting_unit"] == "mapped_identifier"
    # A representative-only FASTA remains accepted without changing its IDs.
    fasta.write_text(">t2\nATGATG\n>t3\nATGATG\n")
    representatives = prepare_genome("z", fasta, gff, tmp_path, **overrides)
    assert [g.gene_id for g in representatives.mapping.genes] == ["t2", "t3"]


def test_ambiguous_locus_is_rejected(tmp_path):
    fasta = write(tmp_path / "cds.fa", ">t1\nATG\n")
    gff = write(tmp_path / "genes.gff", "chr1\tt\tmRNA\t1\t3\t.\t+\t.\tID=t1;Parent=g1,g2\n")
    with pytest.raises(ValueError, match="multiple gene loci"):
        prepare_genome("x", fasta, gff, tmp_path)


@pytest.mark.parametrize(
    "seqid,strand,message",
    [
        ("", "+", "Empty GFF sequence identifier"),
        ("chr 1", "+", "Invalid GFF sequence identifier"),
        ("chr\x00", "+", "Invalid GFF sequence identifier"),
        ("chr1", "INVALID", "Invalid GFF strand"),
    ],
)
def test_invalid_gff_fields_report_the_source_line(tmp_path, seqid, strand, message):
    gff = write(
        tmp_path / "invalid.gff",
        f"##gff-version 3\n\n{seqid}\ttest\tmRNA\t1\t3\t.\t{strand}\t.\tID=t1\n",
    )
    with pytest.raises(ValueError, match=message) as caught:
        annotation_to_genes(gff, {"t1"})
    assert f"{gff}:3" in str(caught.value)


@pytest.mark.parametrize("strands", [("+", "-"), ("+", ".", "-")])
def test_conflicting_gff_segments_are_not_silently_merged(tmp_path, strands):
    gff = write(
        tmp_path / "conflict.gff",
        "".join(
            f"chr1\ttest\tCDS\t{i * 10 + 1}\t{i * 10 + 3}\t.\t{strand}\t0\tParent=t1\n"
            for i, strand in enumerate(strands)
        ),
    )
    with pytest.raises(ValueError, match="conflicting strands") as caught:
        annotation_to_genes(gff, {"t1"})
    assert f"{gff}:{len(strands)}" in str(caught.value)


@pytest.mark.parametrize(
    "strands,expected",
    [
        (("?", "."), "."),
        (("+", "?", "."), "+"),
        (("?", "-"), "-"),
    ],
)
def test_valid_gff_segments_produce_readable_bed(tmp_path, strands, expected):
    fasta = write(tmp_path / "cds.fa", ">t1\nATGATG\n")
    gff = write(
        tmp_path / "genes.gff",
        "".join(
            f"chr1\ttest\tCDS\t{i * 10 + 1}\t{i * 10 + 3}\t.\t{strand}\t0\tParent=t1\n"
            for i, strand in enumerate(strands)
        ),
    )
    prepared = prepare_genome("valid", fasta, gff, tmp_path)
    assert prepared.mapping.feature == "CDS"
    assert prepared.mapping.attribute == "Parent"
    assert prepared.mapping.unresolved_locus_count == 1
    assert read_bed(prepared.bed_path) == prepared.mapping.genes
    (gene,) = prepared.mapping.genes
    assert (gene.start, gene.end, gene.strand) == (0, (len(strands) - 1) * 10 + 3, expected)


@pytest.mark.parametrize(
    "record,message",
    [
        ("chr1\t-1\t5\tt1\n", "Invalid BED interval"),
        ("chr1\t5\t5\tt1\n", "Invalid BED interval"),
        ("\t0\t5\tt1\n", "Empty BED sequence identifier"),
        ("chr1\t0\t5\t\n", "Empty BED gene identifier"),
        ("chr 1\t0\t3\tt1\n", "Invalid BED sequence identifier"),
        ("chr1\t0\t3\tt1\t0\tINVALID\n", "Invalid BED strand"),
    ],
)
def test_invalid_bed_fields_report_the_source_line(tmp_path, record, message):
    bed = write(tmp_path / "invalid.bed", record)
    with pytest.raises(ValueError, match=message) as caught:
        read_bed(bed)
    assert f"{bed}:1" in str(caught.value)

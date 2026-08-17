from pathlib import Path

import pytest

from kffractbias.io import (
    annotation_to_genes,
    parse_attributes,
    read_fasta_ids,
    read_synmap_pairs,
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_parse_attributes_supports_gff3_and_gtf():
    assert parse_attributes("ID=tx1;Name=alpha%20one;Alias=a,b") == {
        "ID": ("tx1",),
        "Name": ("alpha one",),
        "Alias": ("a", "b"),
    }
    assert parse_attributes('gene_id "g1"; transcript_id "t1";') == {
        "gene_id": ("g1",),
        "transcript_id": ("t1",),
    }


def test_annotation_mapping_detects_feature_and_attribute(tmp_path):
    fasta = write(tmp_path / "genes.fa", ">tx1\nATG\n>tx2\nATG\n")
    gff = write(
        tmp_path / "genes.gff3",
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1\t12\t.\t+\t.\tID=tx1;Parent=g1\n"
        "chr1\ttest\tmRNA\t20\t30\t.\t-\t.\tID=tx2;Parent=g1\n",
    )
    mapping = annotation_to_genes(gff, read_fasta_ids(fasta))
    assert mapping.feature == "mRNA"
    assert mapping.attribute == "ID"
    assert mapping.matched_gene_count == 2
    assert [(gene.gene_id, gene.start, gene.end) for gene in mapping.genes] == [
        ("tx1", 0, 12),
        ("tx2", 19, 30),
    ]


def test_annotation_mapping_merges_cds_segments(tmp_path):
    fasta = write(tmp_path / "genes.fa", ">tx1\nATG\n")
    gff = write(
        tmp_path / "genes.gff3",
        "chr1\ttest\tCDS\t1\t3\t.\t+\t0\tParent=tx1\n"
        "chr1\ttest\tCDS\t10\t12\t.\t+\t0\tParent=tx1\n",
    )
    mapping = annotation_to_genes(gff, read_fasta_ids(fasta))
    assert mapping.feature == "CDS"
    assert mapping.attribute == "Parent"
    assert mapping.genes[0].start == 0
    assert mapping.genes[0].end == 12


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


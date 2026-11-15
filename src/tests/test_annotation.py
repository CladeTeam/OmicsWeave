# ruff: noqa: E501
import gzip
import pickle

import numpy as np
import pytest

from omicsweave import errors
from omicsweave.annotation import FEATURE_KINDS, Annotation, IntervalSet
from omicsweave.types import Region, Spliced

GTF = """##format: gtf
chr1\tHAVANA\tgene\t101\t1000\t.\t+\t.\tgene_id "G1"; gene_type "protein_coding"; gene_name "ALPHA"; level 2;
chr1\tHAVANA\ttranscript\t101\t1000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; transcript_type "protein_coding"; transcript_name "ALPHA-201"; tag "basic"; tag "MANE_Select";
chr1\tHAVANA\texon\t101\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number 1;
chr1\tHAVANA\texon\t301\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number 2;
chr1\tHAVANA\tCDS\t321\t400\t.\t+\t0\tgene_id "G1"; transcript_id "T1"; exon_number 2;
chr1\tHAVANA\texon\t901\t1000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number 3;
chr1\tHAVANA\ttranscript\t151\t1000\t.\t+\t.\tgene_id "G1"; transcript_id "T2"; transcript_type "retained_intron"; tag "basic";
chr1\tHAVANA\texon\t151\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T2"; exon_number 1;
chr1\tHAVANA\texon\t901\t1000\t.\t+\t.\tgene_id "G1"; transcript_id "T2"; exon_number 2;
chr1\tHAVANA\tgene\t2001\t3000\t.\t-\t.\tgene_id "G2"; gene_type "lncRNA"; gene_name "BETA"; level 2;
chr1\tHAVANA\texon\t2801\t3000\t.\t-\t.\tgene_id "G2"; transcript_id "T3"; transcript_type "lncRNA";
chr1\tHAVANA\texon\t2001\t2100\t.\t-\t.\tgene_id "G2"; transcript_id "T3"; transcript_type "lncRNA";
chrM\tENSEMBL\texon\t11\t60\t.\t+\t.\tgene_id "G3"; transcript_id "T4"; gene_name "MT-X";
"""


@pytest.fixture
def gtf(tmp_path):
    path = tmp_path / "a.gtf"
    path.write_text(GTF)
    return path


def test_tables(gtf):
    a = Annotation.read(gtf)
    assert a.dialect == "gtf"
    assert a.chroms == ("chr1", "chrM")
    g = a.genes
    assert len(g) == 3
    assert g.id.to_list() == ["G1", "G2", "G3"]
    assert g.name.to_list() == ["ALPHA", "BETA", "MT-X"]
    assert g.biotype.to_list() == ["protein_coding", "lncRNA", ""]
    assert g.start.tolist() == [100, 2000, 10] and g.end.tolist() == [1000, 3000, 60]
    assert g.strand.tolist() == [1, -1, 1]
    assert g.synthesized.tolist() == [False, False, True]
    t = a.transcripts
    assert t.id.to_list() == ["T1", "T2", "T3", "T4"]
    assert t.gene.tolist() == [0, 0, 1, 2]
    assert t.length.tolist() == [300, 350, 300, 50]
    assert t.cds_length.tolist() == [80, 0, 0, 0]
    assert t.has_tag("MANE_Select").tolist() == [True, False, False, False]
    assert t.tags(0) == ["basic", "MANE_Select"]
    assert t.synthesized.tolist() == [False, False, True, True]
    f = a.features
    assert [FEATURE_KINDS[k] for k in f.kind[:4]] == ["exon", "exon", "CDS", "exon"]
    assert len(a.exons) == 8
    assert (a.exons.kind == 0).all()
    assert a.report.is_clean
    assert a.report.synthesized_genes == 1
    assert a.heap_bytes > 0


def test_lookups_and_queries(gtf):
    a = Annotation.read(gtf)
    gene = a.gene("ALPHA")
    assert gene.id == "G1" and gene.region == Region("chr1", 100, 1000)
    assert gene.transcripts.tolist() == [0, 1]
    assert ("level", "2") in gene.attributes
    tx = a.transcript("T1")
    assert tx.tags == ("basic", "MANE_Select") and len(tx.features) == 4
    with pytest.raises(KeyError):
        a.gene("nope")
    assert a.overlaps(("chr1", 950, 2050)).tolist() == [0, 1]
    assert a.overlaps(("chr1", 250, 260), "transcript").tolist() == [0, 1]
    exons = a.overlaps(("chr1", 350, 360), "exon")
    assert a.exons.start[exons].tolist() == [150, 300]
    assert a.overlaps(("chrX", 0, 100)).tolist() == []
    off, rows = a.overlaps_batch([("chr1", 0, 150), ("chrM", 0, 20)], "gene")
    assert off.tolist() == [0, 1, 2] and rows.tolist() == [0, 2]
    rows, d = a.nearest(("chr1", 1500, 1600))
    assert rows.tolist() == [1] and d == 400
    assert a.nearest(("chrX", 1, 2))[1] == -1


def test_subsets_and_views(gtf):
    a = Annotation.read(gtf)
    s = a.subset(transcript_biotypes=["protein_coding"])
    assert s.genes.tolist() == [0] and s.transcripts.tolist() == [0]
    assert a.subset(tags=["basic"]).transcripts.tolist() == [0, 1]
    assert a.subset(genes=["BETA"], chroms=["chr1"]).genes.tolist() == [1]
    with pytest.raises(ValueError):
        a.subset(tags=["no_such_tag"])

    introns = a.introns()
    assert list(zip(introns.start.tolist(), introns.end.tolist(), strict=True)) == [
        (200, 300),
        (400, 900),
        (400, 900),
        (2100, 2800),
    ]
    assert introns.owner.tolist() == [0, 0, 1, 2]
    assert introns.owner_level == "transcript"
    assert len(a.introns(s)) == 2

    exonic = a.exonic_regions()
    assert exonic.start.tolist()[:3] == [100, 900, 2000]
    j = a.junctions()
    assert len(j) == 3
    i = j.intervals.start.tolist().index(400)
    assert j.transcripts[j.transcript_offsets[i] : j.transcript_offsets[i + 1]].tolist() == [0, 1]
    minus = j.intervals.strand.tolist().index(-1)
    assert (j.donor[minus], j.acceptor[minus]) == (2799, 2100)

    bins = a.exon_bins()
    first = [
        (int(s), int(e))
        for s, e in zip(bins.intervals.start[:3], bins.intervals.end[:3], strict=True)
    ]
    assert first == [(100, 150), (150, 200), (200, 300)]
    assert bins.part[:3].tolist() == [1, 2, 3]
    assert bins.transcripts[bins.transcript_offsets[1] : bins.transcript_offsets[2]].tolist() == [
        0,
        1,
    ]

    p = a.promoters(100, 10, level="gene", chrom_lengths={"chr1": 2950, "chrM": 16569})
    assert list(zip(p.start.tolist(), p.end.tolist(), strict=True)) == [
        (0, 110),
        (2950, 2950),
        (0, 20),
    ]
    assert p.name.to_list() == ["G1", "G2", "G3"]
    body = a.gene_bodies(upstream=500)
    assert (body.start[1], body.end[1]) == (2000, 3500)
    tss = a.tss("gene")
    assert tss.start.tolist() == [100, 2999, 10]
    tes = a.tes_windows(10, 20, level="gene")
    assert (tes.start[1], tes.end[1]) == (1980, 2010)
    assert p.overlaps(("chr1", 50, 60)).tolist() == [0]


def test_round_trip_and_bed12(gtf, tmp_path):
    a = Annotation.read(gtf)
    a.write_gtf(tmp_path / "b.gtf")
    assert (tmp_path / "b.gtf").read_text() == GTF
    a.write_gtf(tmp_path / "m.gtf", layout="model")
    m = Annotation.read(tmp_path / "m.gtf")
    assert m.genes.id.to_list() == a.genes.id.to_list()
    assert m.features.start.tolist() == a.features.start.tolist()
    a.write_bed12(tmp_path / "t.bed", a.subset(genes=["G1"]))
    lines = (tmp_path / "t.bed").read_text().splitlines()
    assert lines[0] == "chr1\t100\t1000\tT1\t0\t+\t320\t400\t0\t3\t100,100,100,\t0,200,800,"
    none = Annotation.read(gtf, attributes=None)
    with pytest.raises(errors.UnsupportedError):
        none.write_gtf(tmp_path / "x.gtf", layout="original")


def test_options_compression_and_pickle(gtf, tmp_path):
    gz = tmp_path / "a.gtf.gz"
    gz.write_bytes(gzip.compress(GTF.encode()))
    a = Annotation.read(gz, aliases={"chrM": "MT"}, threads=2)
    assert a.chroms == ("chr1", "MT")
    assert a.match_chroms(["1", "chrM", "chr2"]).tolist() == [0, 1, -1]
    only = Annotation.read(gtf, attributes=["level"])
    assert only.gene("G1").attributes == (("level", "2"),)
    assert only.transcript("T1").attributes == ()
    b = pickle.loads(pickle.dumps(a))
    assert b.chroms == a.chroms and b.genes.id.to_list() == a.genes.id.to_list()


def test_problems(tmp_path):
    bad = tmp_path / "bad.gtf"
    bad.write_text(
        'c\ts\ttranscript\t1\t50\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
        'c\ts\texon\t1\t60\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
    )
    a = Annotation.read(bad)
    assert a.report.count("outside_parent") == 1
    assert a.report.problems[0].line == 2 and a.report.problems[0].id == "t"
    with pytest.raises(errors.FormatError) as e:
        Annotation.read(bad, strict=True)
    assert e.value.line == 2
    broken = tmp_path / "broken.gtf"
    broken.write_text('c\ts\texon\t0\t9\t.\t+\t.\tgene_id "g";\n')
    with pytest.raises(errors.FormatError):
        Annotation.read(broken)


def test_gff3(tmp_path):
    path = tmp_path / "a.gff3"
    path.write_text(
        "##gff-version 3\n"
        "c\tsrc\tgene\t1\t500\t.\t-\t.\tID=gene-A;Name=A;gene_biotype=protein_coding\n"
        "c\tsrc\tmRNA\t1\t500\t.\t-\t.\tID=rna-A1;Parent=gene-A;tag=MANE Select\n"
        "c\tsrc\texon\t1\t100\t.\t-\t.\tParent=rna-A1\n"
        "c\tsrc\texon\t401\t500\t.\t-\t.\tParent=rna-A1\n"
    )
    a = Annotation.read(path)
    assert a.dialect == "gff3"
    assert a.genes.name.to_list() == ["A"]
    assert a.transcripts.has_tag("MANE_Select").tolist() == [True]
    assert a.introns().start.tolist() == [100]


def test_interval_sets_from_arrays():
    s = IntervalSet.from_arrays(
        ["chr1", "chr2"],
        np.array([0, 1, 0]),
        np.array([10, 0, 50]),
        np.array([20, 5, 60]),
        names=["p1", "p2", "p3"],
    )
    assert len(s) == 3 and s.chroms == ("chr1", "chr2")
    assert s.overlaps(("chr1", 15, 55)).tolist() == [0, 2]
    assert s.name.to_list() == ["p1", "p2", "p3"]
    assert s.owner_level is None
    off, rows = s.overlaps_batch([("chr2", 0, 1), ("chr3", 0, 1)])
    assert off.tolist() == [0, 1, 1] and rows.tolist() == [1]
    with pytest.raises(errors.InvalidRegionError):
        IntervalSet.from_arrays(["c"], np.array([0]), np.array([5]), np.array([4]))


MINUS = """chr2\tX\texon\t5001\t5100\t.\t-\t.\tgene_id "G5"; transcript_id "T5";
chr2\tX\texon\t5201\t5300\t.\t-\t.\tgene_id "G5"; transcript_id "T5";
chr2\tX\tCDS\t5221\t5300\t.\t-\t0\tgene_id "G5"; transcript_id "T5";
chr2\tX\tstop_codon\t5218\t5220\t.\t-\t0\tgene_id "G5"; transcript_id "T5";
"""


def test_transcript_parts_as_spliced_keys(gtf, tmp_path):
    a = Annotation.read(gtf)
    t1 = a.locus("T1")
    assert t1 == Spliced("chr1", [(100, 200), (300, 400), (900, 1000)], "+", "T1")
    assert len(t1) == a.transcripts.length[0]
    assert a.locus("T1", "cds").spans == ((320, 400),)
    assert a.locus("T1", "utr5").spans == ((100, 200), (300, 320))
    assert a.locus("T1", "utr3").spans == ((900, 1000),)
    assert a.locus(2).strand == "-" and a.locus(2, "cds").spans == ()
    assert [k.name for k in a.loci()] == ["T1", "T2", "T3", "T4"]
    assert [len(k) for k in a.loci([0, 3], "exons")] == [300, 50]
    path = tmp_path / "m.gtf"
    path.write_text(MINUS)
    m = Annotation.read(path)
    # The stop codon, listed apart from CDS, belongs to the coding part.
    assert m.locus("T5", "cds").spans == ((5217, 5300),)
    assert m.locus("T5", "utr5").spans == ()
    assert m.locus("T5", "utr3").spans == ((5000, 5100), (5200, 5217))
    with pytest.raises(ValueError):
        a.locus("T1", "intron")

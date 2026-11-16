"""Options the Rust layer has that every reader, writer and conversion
exposes: writer threads and block layouts, sidecars read back, 2bit
writing, CRAM decoding, cache keys, stranded caches and contact maps."""

import numpy as np
import pytest

import omicsweave as ow
from omicsweave import errors
from omicsweave.annotation import Annotation, read_aliases
from omicsweave.cache import BlockCache
from omicsweave.formats import (
    BbiSidecar,
    BigWigWriter,
    HicBlocks,
    HicBuilder,
    HicHeader,
    HicSidecar,
    TwoBitWriter,
    read_sidecar,
)
from omicsweave.query import (
    AlignmentReader,
    ContactReader,
    SequenceReader,
    StrandedTrackReader,
    TrackReader,
)
from omicsweave.types import ContactWindow

CHROMS = {"chr1": 50_000, "chr2": 23_000}


def track(path, rng, **kw):
    starts = np.arange(0, 40_000, 7, dtype=np.uint32)
    with BigWigWriter(path, CHROMS, **kw) as w:
        w.add("chr1", starts, starts + 5, rng.normal(size=len(starts)).astype(np.float32))
    return path


def contacts(path, rng, **kw):
    n = 4000
    c1, c2 = rng.integers(0, 2, n), rng.integers(0, 2, n)
    lens = np.array(list(CHROMS.values()))
    header = HicHeader(genome="g", chroms=CHROMS, bp_resolutions=(1000, 5000))
    with HicBuilder(path, header, **kw) as b:
        b.add(c1, rng.integers(0, lens[c1]), c2, rng.integers(0, lens[c2]), np.ones(n, np.float32))
    return path


def test_writer_threads_do_not_change_the_file(tmp_path):
    one = track(tmp_path / "one.bw", np.random.default_rng(1))
    four = track(tmp_path / "four.bw", np.random.default_rng(1), threads=4)
    assert one.read_bytes() == four.read_bytes()
    a = contacts(tmp_path / "one.hic", np.random.default_rng(2))
    b = contacts(tmp_path / "four.hic", np.random.default_rng(2), threads=4)
    assert a.read_bytes() == b.read_bytes()
    with pytest.raises(ValueError, match="at least 1"):
        track(tmp_path / "zero.bw", np.random.default_rng(1), threads=0)


def test_hic_block_layouts_read_the_same(tmp_path):
    default = contacts(tmp_path / "d.hic", np.random.default_rng(3))
    juicer = contacts(
        tmp_path / "j.hic", np.random.default_rng(3), blocks=HicBlocks("juicer", "juicer")
    )
    fixed = contacts(tmp_path / "f.hic", np.random.default_rng(3), blocks=HicBlocks(16, 16))
    assert len({default.read_bytes(), juicer.read_bytes(), fixed.read_bytes()}) == 3
    w = ContactWindow(("chr1", 0, 50_000), ("chr2", 0, 23_000), 1000)
    with ContactReader(default) as d, ContactReader(juicer) as j, ContactReader(fixed) as f:
        np.testing.assert_array_equal(d.matrix(w).values, j.matrix(w).values)
        np.testing.assert_array_equal(d.matrix(w).values, f.matrix(w).values)
        assert d.maps == (("chr1", "chr1"), ("chr1", "chr2"), ("chr2", "chr2"))


def test_sidecars_read_back(tmp_path):
    bw = track(tmp_path / "t.bw", np.random.default_rng(4), sidecar=True)
    s = read_sidecar(bw)
    assert isinstance(s, BbiSidecar) and s.levels[0].reduction == 0
    base = s.levels[0].blocks
    assert len(s.levels[0]) > 0 and (base["inflated"] >= base["size"]).all()
    assert all(level.reduction > 0 for level in s.levels[1:])
    hic = contacts(tmp_path / "c.hic", np.random.default_rng(5), sidecar=True)
    h = read_sidecar(hic)
    assert isinstance(h, HicSidecar)
    assert {(m.chrom1, m.chrom2, m.unit, m.resolution) for m in h.maps} >= {(1, 1, "BP", 1000)}
    assert sum(int(m.blocks["records"].sum()) for m in h.maps if m.resolution == 1000) > 0
    with pytest.raises(errors.OmicsWeaveError):
        read_sidecar(track(tmp_path / "none.bw", np.random.default_rng(4)))


def test_two_bit_writer_round_trips(tmp_path):
    seqs = {"a": "ACGTNNNNacgtAC", "b": "GGGGTTTTnnA"}
    with TwoBitWriter(tmp_path / "g.2bit", seqs) as w:
        for name, seq in seqs.items():
            assert w.next_name == name
            w.add(name, seq)
        assert w.next_name is None
    with SequenceReader(tmp_path / "g.2bit", soft_mask=True) as r:
        assert dict(r.chroms) == {k: len(v) for k, v in seqs.items()}
        assert r.fetch(("a", 0, 14)).text == "ACGTNNNNacgtAC"
        assert r.fetch(("b", 0, 11)).text == "GGGGTTTTNNA"
    w = TwoBitWriter(tmp_path / "short.2bit", ["a", "b"])
    w.add("a", "ACGT")
    with pytest.raises(errors.OmicsWeaveError):
        w.close()


def test_stranded_tracks_share_a_cache_and_clone(tmp_path):
    plus = track(tmp_path / "p.bw", np.random.default_rng(6))
    minus = track(tmp_path / "m.bw", np.random.default_rng(7))
    cache = BlockCache(max_bytes=1 << 24)
    with StrandedTrackReader(plus, minus, cache=cache) as rna:
        rna.signal(("chr1", 0, 10_000, "-"))
        loads = cache.stats()["loads"]
        assert loads > 0
        clone = rna.clone()
        np.testing.assert_array_equal(
            clone.signal(("chr1", 0, 10_000, "-")).values,
            rna.signal(("chr1", 0, 10_000, "-")).values,
        )
        assert cache.stats()["loads"] == loads
        clone.close()
        assert not rna.closed
        with TrackReader(plus) as t:
            assert rna.plus.summary == t.summary


def test_cram_decoding_options_belong_to_cram(tmp_path):
    sam = tmp_path / "r.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:100\nr\t0\tc\t5\t9\t4M\t*\t0\t0\tACGT\t*\n"
    )
    bam = ow.convert.run(sam, tmp_path / "r.bam", index=True, threads=1).target
    with AlignmentReader(bam) as reads:
        assert len(reads.fetch(("c", 0, 100))) == 1
    with pytest.raises(errors.InvalidError, match="CRAM decoding"):
        AlignmentReader(bam, verify_checksums=False)


def test_conversion_options_are_checked_when_planned(tmp_path):
    (tmp_path / "a.fa").write_bytes(b">a\nACGT\n")
    with pytest.raises(errors.InvalidError, match="index"):
        ow.convert.plan(tmp_path / "a.fa", tmp_path / "a.2bit", index=True)
    with pytest.raises(errors.InvalidError, match="level"):
        ow.convert.plan(tmp_path / "a.fa", tmp_path / "a.2bit", level=3)
    with pytest.raises(errors.InvalidError, match="verify_checksums"):
        ow.convert.plan(tmp_path / "a.fa", tmp_path / "a.fa.gz", verify_checksums=False)
    sam = tmp_path / "r.sam"
    sam.write_text("@SQ\tSN:c\tLN:100\nr\t0\tc\t5\t9\t4M\t*\t0\t0\tACGT\t*\n")
    with pytest.raises(errors.InvalidError, match="min_shift"):
        ow.convert.plan(sam, tmp_path / "r.bam", index="bai", min_shift=12)
    report = ow.convert.run(sam, tmp_path / "r.bam", index="csi", min_shift=12, threads=1)
    assert [f.format for f in report.companions] == ["csi"]
    (tmp_path / "a.bedgraph").write_text("c\t0\t10\t1.5\nc\t10\t20\t2\n")
    small = ow.convert.run(
        tmp_path / "a.bedgraph", tmp_path / "s.bw", bbi=ow.formats.BbiOptions(zoom_levels=0)
    )
    with TrackReader(small.target) as t:
        assert t.zoom_reductions == ()


def test_annotation_others_containment_and_aliases(tmp_path):
    gff = tmp_path / "a.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tsrc\tregion\t1\t1000\t.\t+\t.\tID=chr1\n"
        "chr1\tsrc\tgene\t101\t500\t.\t+\t.\tID=g1;gene_name=G1\n"
        "chr1\tsrc\tmRNA\t101\t500\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tsrc\texon\t101\t200\t.\t+\t.\tParent=t1\n"
        "chr1\tsrc\texon\t301\t500\t.\t+\t.\tParent=t1\n"
    )
    ann = Annotation.read(gff, threads=1)
    assert ann.others.type.to_list() == ["region"] and ann.others.id.to_list() == ["chr1"]
    exons = ann.exonic_regions()
    assert exons.overlaps(("chr1", 150, 350)).tolist() == [0, 1]
    assert exons.overlaps(("chr1", 50, 250), mode="within").tolist() == [0]
    assert exons.overlaps(("chr1", 120, 130), mode="containing").tolist() == [0]
    (tmp_path / "aliases").write_text("# from to\n1 chr1\n\nMT chrM\n")
    assert read_aliases(tmp_path / "aliases") == {"1": "chr1", "MT": "chrM"}
    (tmp_path / "bad").write_text("1\n")
    with pytest.raises(errors.FormatError):
        read_aliases(tmp_path / "bad")

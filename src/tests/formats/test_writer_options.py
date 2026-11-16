import numpy as np
import pytest

from omicsweave import errors
from omicsweave.formats import (
    BbiSidecar,
    HicBlocks,
    HicBuilder,
    HicHeader,
    HicSidecar,
    TwoBitWriter,
    read_sidecar,
)
from omicsweave.query import (
    ContactReader,
    SequenceReader,
)
from omicsweave.types import ContactWindow

from ..support.tracks import CHROMS, track


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

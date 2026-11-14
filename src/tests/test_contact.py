import numpy as np
import pytest

from omicsweave import errors
from omicsweave.formats import Hic, HicBuilder, HicHeader, HicWriter, Sidecar
from omicsweave.query import ContactReader
from omicsweave.types import ContactWindow, Region

CHROMS = {"chr1": 50_000, "chr2": 23_000}
RESOLUTIONS = (1000, 5000)


def random_contacts(rng, n=3000):
    names = list(CHROMS)
    c1 = rng.integers(0, 2, n)
    c2 = rng.integers(0, 2, n)
    lens = np.array([CHROMS[k] for k in names])
    p1 = rng.integers(0, lens[c1])
    p2 = rng.integers(0, lens[c2])
    count = rng.integers(1, 4, n).astype(np.float32)
    return c1, p1, c2, p2, count


def expected(contacts, q: ContactWindow):
    names = list(CHROMS)
    rows, cols = q.row_bins, q.col_bins
    out = np.zeros(q.shape, dtype=np.float64)
    for c1, p1, c2, p2, n in zip(*contacts, strict=True):
        a, b = names[c1], names[c2]
        x, y = int(p1) // q.resolution, int(p2) // q.resolution
        cells = {(a, x, b, y), (b, y, a, x)}
        for ra, rx, cb, cy in cells:
            if ra == q.rows.chrom and cb == q.cols.chrom and rx in rows and cy in cols:
                out[rx - rows.start, cy - cols.start] += n
    return out


def header():
    return HicHeader(genome="test", chroms=CHROMS, bp_resolutions=RESOLUTIONS)


def build(path, contacts, **kw):
    with HicBuilder(path, header(), **kw) as b:
        b.add(*contacts)


def random_window(rng):
    names = list(CHROMS)
    res = RESOLUTIONS[int(rng.integers(len(RESOLUTIONS)))]

    def region():
        name = names[int(rng.integers(2))]
        start = int(rng.integers(0, CHROMS[name]))
        return Region(name, start, int(rng.integers(start, CHROMS[name] + 1)))

    return ContactWindow(region(), region(), res)


def test_built_map_matches_binned_contacts(tmp_path):
    rng = np.random.default_rng(0)
    contacts = random_contacts(rng)
    path = tmp_path / "m.hic"
    build(path, contacts)
    assert not (tmp_path / "m.hic.sidecar").exists()
    f = ContactReader(path)
    assert f.genome == "test" and f.version == 9
    assert {k: v for k, v in f.chroms.items() if k in CHROMS} == CHROMS
    assert f.resolutions() == RESOLUTIONS
    assert f.resolutions("FRAG") == ()
    for _ in range(25):
        q = random_window(rng)
        m = f.matrix(q)
        assert m.values.dtype == np.float32 and m.values.shape == q.shape
        np.testing.assert_allclose(m.values, expected(contacts, q))
        cells = f.pixels(q)
        dense = np.zeros(q.shape, dtype=np.float32)
        rows = cells.row // q.resolution - q.row_bins.start
        cols = cells.col // q.resolution - q.col_bins.start
        dense[rows, cols] = cells.value
        np.testing.assert_array_equal(dense, m.values)
        assert np.all(cells.value != 0)


def test_writer_round_trip(tmp_path):
    rng = np.random.default_rng(1)
    path = tmp_path / "w.hic"
    x = rng.integers(0, 50, 200).astype(np.uint32)
    y = rng.integers(0, 50, 200).astype(np.uint32)
    v = rng.random(200).astype(np.float32) + 0.5
    with HicWriter(path, header()) as w:
        w.write_matrix("chr1", "chr1", 1000, x, y, v)
    f = ContactReader(path)
    q = ContactWindow(("chr1", 0, 50_000), ("chr1", 0, 50_000), 1000)
    want = np.zeros((50, 50))
    for a, b, n in zip(x, y, v, strict=True):
        lo, hi = min(a, b), max(a, b)
        want[lo, hi] += n
    want = want + np.triu(want, 1).T
    np.testing.assert_allclose(f.matrix(q).values, want, rtol=1e-6)


def test_out_and_errors(tmp_path):
    path = tmp_path / "m.hic"
    build(path, random_contacts(np.random.default_rng(2), 100))
    f = ContactReader(path)
    q = ContactWindow(("chr1", 0, 10_000), ("chr2", 0, 3000), 1000)
    out = np.empty(q.shape, dtype=np.float32)
    assert f.matrix(q, out=out).values is out
    with pytest.raises(errors.UnknownChromError):
        f.matrix(ContactWindow(("chrX", 0, 10), ("chr1", 0, 10), 1000))
    with pytest.raises(errors.InvalidError):
        ContactWindow(("chr1", 0, 10), ("chr1", 0, 10), 0)
    with pytest.raises(errors.OmicsWeaveError):
        f.matrix(ContactWindow(("chr1", 0, 10), ("chr1", 0, 10), 777))


def test_sidecar_written_equals_built(tmp_path):
    contacts = random_contacts(np.random.default_rng(3), 500)
    path, sidecar = tmp_path / "m.hic", tmp_path / "m.hic.sidecar"
    build(path, contacts, sidecar=True)
    written = sidecar.read_bytes()
    sidecar.unlink()
    assert Hic(path).build_sidecar() == Sidecar(sidecar)
    assert sidecar.read_bytes() == written

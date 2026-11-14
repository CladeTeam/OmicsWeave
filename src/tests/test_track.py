import numpy as np
import pytest

from omicsweave import errors
from omicsweave.formats import BbiOptions, BigWig, BigWigSummary, BigWigWriter, Sidecar
from omicsweave.query import IntervalReader, TrackReader

CHROMS = {"chr1": 5000, "chr2": 1200, "chrM": 300}


def random_items(rng):
    """Sorted, non-overlapping items per chromosome."""
    items = {}
    for name, length in CHROMS.items():
        cuts = np.unique(rng.integers(0, length + 1, int(rng.integers(2, 80))))
        starts, ends = cuts[:-1], cuts[1:]
        keep = rng.random(len(starts)) < 0.7
        values = rng.normal(size=len(starts)).astype(np.float32)
        items[name] = (starts[keep], ends[keep], values[keep])
    return items


def write_track(path, items, **kw):
    with BigWigWriter(path, CHROMS, **kw) as w:
        for name, (s, e, v) in items.items():
            w.add(name, s, e, v)


def dense(items, name, start, end):
    out = np.full(end - start, np.nan, dtype=np.float32)
    for s, e, v in zip(*items[name], strict=True):
        lo, hi = max(int(s), start), min(int(e), end)
        if lo < hi:
            out[lo - start : hi - start] = v
    return out


@pytest.mark.parametrize("seed", range(4))
def test_values_and_intervals_match_written_items(tmp_path, seed):
    rng = np.random.default_rng(seed)
    items = random_items(rng)
    path = tmp_path / "t.bw"
    options = BbiOptions(items_per_slot=int(rng.integers(1, 20)), block_size=4)
    write_track(path, items, options=options)
    f = TrackReader(path)
    assert dict(f.chroms) == CHROMS
    for _ in range(40):
        name = list(CHROMS)[int(rng.integers(len(CHROMS)))]
        start = int(rng.integers(0, CHROMS[name] + 50))
        end = int(rng.integers(start, CHROMS[name] + 100))
        got = f.signal((name, start, end))
        np.testing.assert_array_equal(got.values, dense(items, name, start, end))
        iv = f.intervals((name, start, end))
        rebuilt = np.full(end - start, np.nan, dtype=np.float32)
        for s, e, v in zip(iv.start, iv.end, iv.value, strict=True):
            assert start <= s < e <= end
            rebuilt[s - start : e - start] = v
        np.testing.assert_array_equal(rebuilt, got.values)


def test_summary(tmp_path):
    items = {"chr1": (np.array([0, 10]), np.array([5, 20]), np.array([1.0, 2.0], np.float32))}
    path = tmp_path / "t.bw"
    write_track(path, items)
    s = TrackReader(path).summary
    assert s == BigWigSummary(15, 1.0, 2.0, 25.0, 45.0)


def test_sidecars_are_explicit(tmp_path):
    items = random_items(np.random.default_rng(9))
    path, sidecar = tmp_path / "t.bw", tmp_path / "t.bw.sidecar"
    write_track(path, items)
    assert not sidecar.exists()
    TrackReader(path).signal(("chr1", 0, 100))
    assert not sidecar.exists()
    assert BigWig(path).build_sidecar() == Sidecar(sidecar)
    built = sidecar.read_bytes()
    sidecar.unlink()
    write_track(path, items, sidecar=True)
    assert sidecar.read_bytes() == built


def test_out_argument_and_errors(tmp_path):
    path = tmp_path / "t.bw"
    write_track(path, {"chr1": (np.array([0]), np.array([10]), np.array([3.0], np.float32))})
    f = TrackReader(path)
    out = np.zeros(20, dtype=np.float32)
    assert f.signal(("chr1", 0, 20), out=out).values is out
    assert np.isnan(out[10:]).all() and (out[:10] == 3).all()
    with pytest.raises(errors.UnknownChromError):
        f.signal(("chrX", 0, 1))
    with pytest.raises(errors.InvalidError):
        IntervalReader(path)
    with BigWigWriter(tmp_path / "u.bw", CHROMS) as w:
        with pytest.raises(errors.UnknownChromError):
            w.add("chrX", 0, 1, 1.0)
        w.add("chr1", 10, 20, 1.0)
        with pytest.raises(errors.InvalidError):
            w.add("chr1", 5, 8, 1.0)

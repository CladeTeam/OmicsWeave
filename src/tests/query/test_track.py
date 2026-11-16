import pickle

import numpy as np
import pytest

from omicsweave import errors
from omicsweave.formats import BbiOptions, BigWig, BigWigSummary, BigWigWriter, Sidecar
from omicsweave.query import IntervalReader, StrandedTrackReader, TrackReader
from omicsweave.types import Region, Spliced, Stat

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


def brute_bins(values, resolution, phase, stat):
    """Per-base reference of one bin row (values already 5' to 3')."""
    n = 0 if len(values) == 0 else -(-(len(values) + phase) // resolution)
    out = np.empty(n)
    for k in range(n):
        lo, hi = max(k * resolution - phase, 0), min((k + 1) * resolution - phase, len(values))
        part = values[lo:hi].astype(np.float64)
        data = part[~np.isnan(part)]
        if stat == "coverage":
            out[k] = len(data) / len(part)
        elif len(data) == 0:
            out[k] = np.nan
        elif stat == "mean":
            out[k] = data.mean()
        elif stat == "sum":
            out[k] = data.sum()
        elif stat == "min":
            out[k] = data.min()
        elif stat == "max":
            out[k] = data.max()
        else:
            out[k] = 0.0 if len(data) == 1 else data.std(ddof=1)
    return out


@pytest.mark.parametrize("seed", range(4))
def test_bins_summarize_signal_on_any_key(tmp_path, seed):
    rng = np.random.default_rng(seed)
    items = random_items(rng)
    path = tmp_path / "t.bw"
    write_track(path, items)
    f = TrackReader(path)
    stats = ("mean", "min", "max", "sum", "std", "coverage")
    for _ in range(30):
        name = list(CHROMS)[int(rng.integers(len(CHROMS)))]
        strand = str(rng.choice(["+", "-", "."]))
        if rng.random() < 0.5:
            start = int(rng.integers(0, CHROMS[name]))
            key = Region(name, start, start + int(rng.integers(0, 900)), strand)
        else:
            cuts = np.sort(rng.integers(0, CHROMS[name], 2 * int(rng.integers(1, 4))))
            key = Spliced(name, cuts.reshape(-1, 2), strand)
        resolution = int(rng.integers(1, 200))
        phase = int(rng.integers(0, resolution))
        bins = f.bins(key, resolution, stats, phase=phase)
        values = f.signal(key).values
        assert bins.values.shape == (6, len(bins))
        assert bins.edges[0] == 0 and bins.edges[-1] == len(key)
        for stat in stats:
            np.testing.assert_allclose(
                bins[stat], brute_bins(values, resolution, phase, stat), rtol=1e-5, atol=1e-4
            )


def test_bins_options_and_out(tmp_path):
    path = tmp_path / "t.bw"
    write_track(path, random_items(np.random.default_rng(1)))
    f = TrackReader(path)
    one = f.bins(("chr1", 0, 1000), 100)
    assert one.stats == (Stat.MEAN,) and one.values.shape == (1, 10)
    out = np.zeros((2, 10), dtype=np.float32)
    got = f.bins(("chr1", 0, 1000), 100, [Stat.MAX, "coverage"], out=out)
    assert got.values is out
    np.testing.assert_array_equal(out[0], f.bins(("chr1", 0, 1000), 100, "max").values[0])
    zoomed = f.bins(("chr1", 0, 4000), 1000, ("mean", "coverage"), exact=False)
    assert zoomed.values.shape == (2, 4)
    for bad in (dict(phase=100), dict(phase=-1)):
        with pytest.raises(errors.InvalidError):
            f.bins(("chr1", 0, 1000), 100, **bad)
    with pytest.raises(errors.InvalidError):
        f.bins(("chr1", 0, 1000), 0)
    with pytest.raises(ValueError):
        f.bins(("chr1", 0, 1000), 100, "median")
    with pytest.raises(errors.InvalidError):
        f.bins(("chr1", 0, 1000), 100, ("mean", "mean"))
    with pytest.raises(ValueError):
        f.bins(("chr1", 0, 1000), 100, out=np.zeros((10, 1), dtype=np.float32))
    assert len(f.bins(("chr1", 10, 10), 100)) == 0


def test_reverse_and_spliced_signal(tmp_path):
    rng = np.random.default_rng(5)
    items = random_items(rng)
    path = tmp_path / "t.bw"
    write_track(path, items)
    f = TrackReader(path)
    forward = f.signal(("chr1", 100, 900)).values
    np.testing.assert_array_equal(f.signal(("chr1", 100, 900, "-")).values, forward[::-1])
    key = Spliced("chr1", [(100, 200), (500, 520), (520, 600)], "-", name="t1")
    want = np.concatenate([dense(items, "chr1", 100, 200), dense(items, "chr1", 500, 600)])[::-1]
    np.testing.assert_array_equal(f.signal(key).values, want)
    assert f.signal(key).locus is key


def test_stranded_track_reads_the_key_strand(tmp_path):
    rng = np.random.default_rng(9)
    plus, minus = random_items(rng), random_items(rng)
    write_track(tmp_path / "p.bw", plus)
    write_track(tmp_path / "m.bw", minus)
    with StrandedTrackReader(tmp_path / "p.bw", tmp_path / "m.bw") as f:
        np.testing.assert_array_equal(
            f.signal(("chr1", 0, 500, "+")).values, dense(plus, "chr1", 0, 500)
        )
        np.testing.assert_array_equal(
            f.signal(("chr1", 0, 500, "-")).values, dense(minus, "chr1", 0, 500)[::-1]
        )
        assert f.bins(("chr1", 0, 500, "-"), 50).values.shape == (1, 10)
        with pytest.raises(errors.InvalidError):
            f.signal(("chr1", 0, 500))
    pickle.loads(pickle.dumps(f))

import gzip
from types import SimpleNamespace

import numpy as np
import pytest

import omicsweave as ow
from omicsweave import _native
from omicsweave.formats import BbiOptions, Bed, BigBed, BigBedSchema, BigBedWriter, GzipBed
from omicsweave.query import IntervalReader
from omicsweave.types import StringVector

CHROMS = {"chr1": 4000, "chr2": 900}


def random_features(rng):
    features = []
    for name, length in CHROMS.items():
        n = int(rng.integers(1, 120))
        starts = np.sort(rng.integers(0, length, n))
        ends = np.minimum(starts + rng.integers(1, 400, n), length)
        for s, e in zip(starts, ends, strict=True):
            note = f"f{len(features)}\t{int(rng.integers(0, 1000))}\t{'+-'[int(rng.integers(2))]}"
            features.append((name, int(s), int(e), note))
    return features


def overlapping(features, name, start, end):
    return [
        (s, e, note)
        for c, s, e, note in features
        if c == name and start < end and e > start and s < end
    ]


def check(f, features, rng):
    for _ in range(60):
        name = list(CHROMS)[int(rng.integers(len(CHROMS)))]
        start = int(rng.integers(0, CHROMS[name]))
        end = int(rng.integers(start, CHROMS[name] + 50))
        got = f.fetch((name, start, end))
        assert got.start.dtype == np.uint32 and got.end.dtype == np.uint32
        rows = list(
            zip(got.start.tolist(), got.end.tolist(), got.annotation.to_list(), strict=True)
        )
        assert rows == overlapping(features, name, start, end)


@pytest.mark.parametrize("seed", range(3))
def test_bigbed_round_trip(tmp_path, seed):
    rng = np.random.default_rng(seed)
    features = random_features(rng)
    path = tmp_path / "f.bb"
    schema = BigBedSchema(field_count=6, auto_sql='table x\n"x"\n()\n')
    options = BbiOptions(items_per_slot=int(rng.integers(1, 30)), block_size=3)
    with BigBedWriter(path, CHROMS, schema, options) as w:
        for name in CHROMS:
            rows = [r for r in features if r[0] == name]
            w.add(name, [r[1] for r in rows], [r[2] for r in rows], [r[3] for r in rows])
    f = IntervalReader(path)
    assert f.source == BigBed(path)
    assert dict(f.chroms) == CHROMS
    assert f.schema == BigBedSchema(6, 6, 'table x\n"x"\n()\n')
    check(f, features, rng)


@pytest.mark.parametrize("compressed", [False, True])
def test_bed(tmp_path, compressed):
    rng = np.random.default_rng(7)
    features = random_features(rng)
    order = rng.permutation(len(features))
    text = "".join(
        f"{features[i][0]}\t{features[i][1]}\t{features[i][2]}\t{features[i][3]}\n" for i in order
    )
    path = tmp_path / ("f.bed.gz" if compressed else "f.bed")
    path.write_bytes(gzip.compress(text.encode()) if compressed else text.encode())
    f = ow.open(path)
    assert isinstance(f, IntervalReader)
    assert f.source == (GzipBed(path) if compressed else Bed(path))
    assert f.schema is None
    assert set(f.chroms) == set(CHROMS) and set(f.chroms.values()) == {None}
    stable = sorted(features, key=lambda r: (r[0], r[1]))
    for _ in range(40):
        name = list(CHROMS)[int(rng.integers(len(CHROMS)))]
        start = int(rng.integers(0, CHROMS[name]))
        end = int(rng.integers(start, CHROMS[name] + 50))
        got = f.fetch((name, start, end))
        rows = sorted(
            zip(got.start.tolist(), got.end.tolist(), got.annotation.to_list(), strict=True)
        )
        assert rows == sorted(overlapping(stable, name, start, end))
        assert np.all(np.diff(got.start.astype(np.int64)) >= 0)


def test_string_vector():
    s = StringVector.from_list(["a", "", "héllo"])
    assert len(s) == 3
    assert s[2] == "héllo" and s[-1] == "héllo"
    assert s[0:2] == ["a", ""]
    assert s.bytes(2) == "héllo".encode()
    assert s.offsets.dtype == np.int64 and s.offsets.tolist() == [0, 1, 1, 7]
    with pytest.raises(IndexError):
        s[3]


@pytest.mark.parametrize("n", [0, 1, 257])
@pytest.mark.parametrize("annotation", ["", "héllo"])
def test_scalar_annotation_matches_string_column(tmp_path, n, annotation):
    starts = np.arange(n, dtype=np.uint32) * 2
    ends = starts + 1
    for name, notes in [("scalar", annotation), ("column", [annotation] * n)]:
        path = tmp_path / f"{name}.bb"
        with BigBedWriter(
            path, {"c": 1000}, BigBedSchema(field_count=4), BbiOptions(zoom_levels=0)
        ) as writer:
            writer.add("c", starts, ends, notes)
        with IntervalReader(path) as reader:
            result = reader.fetch(("c", 0, 1000))
            np.testing.assert_array_equal(result.start, starts)
            np.testing.assert_array_equal(result.end, ends)
            np.testing.assert_array_equal(
                result.annotation.offsets,
                np.arange(n + 1, dtype=np.int64) * len(annotation.encode()),
            )
            assert result.annotation.to_list() == [annotation] * n


def test_scalar_annotation_owns_writable_buffers(tmp_path, monkeypatch):
    columns = []
    monkeypatch.setattr(
        _native,
        "BigBedWriter",
        lambda *args: SimpleNamespace(add=lambda *args: columns.append(args)),
    )
    writer = BigBedWriter(tmp_path / "buffers.bb", {"c": 100})
    writer.add("c", [1, 3], [2, 4], "héllo")
    offsets, data = columns[0][-2:]
    assert offsets.flags.owndata and offsets.flags.writeable
    assert data.flags.owndata and data.flags.writeable
    np.testing.assert_array_equal(offsets, [0, 6, 12])
    assert data.tobytes() == "héllohéllo".encode()


def test_empty_scalar_annotation_does_not_encode(tmp_path):
    with BigBedWriter(tmp_path / "empty.bb", {"c": 100}) as writer:
        writer.add("c", [], [], "\ud800")

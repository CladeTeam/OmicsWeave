"""A cache under CRAM queries: every field selection reads the same records."""

import gzip
import multiprocessing
import pickle
import shutil

import numpy as np
import pytest

import omicsweave as ow
from omicsweave.query import AlignmentReader

from ..support.data import CRAM_DATA as DATA


@pytest.fixture(params=["ref", "embed"])
def source(tmp_path, request):
    name = request.param
    for filename in ("ref.fa", f"{name}.cram", f"{name}.cram.crai"):
        shutil.copy(DATA / filename, tmp_path / filename)
    reference = ow.formats.Fasta(tmp_path / "ref.fa")
    reference.build_index()
    return tmp_path / f"{name}.cram", reference if name == "ref" else None


def snapshot(reader, fields="all", region=("chr1", 0, 5000)):
    result = reader.fetch(region, fields)
    out = [result.start.tobytes(), result.end.tobytes(), result.flag.tobytes()]
    for name in ("names", "seqs", "quals", "aux", "cigars"):
        column = getattr(result, name)
        out.append(None if column is None else column.offsets.tobytes())
        if column is not None:
            out.append((column.ops if name == "cigars" else column.data).tobytes())
    return tuple(out)


def test_cram_cache_serves_every_field_selection_and_clones(source):
    path, reference = source
    cache = ow.cache.BlockCache(max_bytes=16 << 20)
    with AlignmentReader(path, reference=reference) as plain:
        expected = {f: snapshot(plain, f) for f in ["all", (), ("seq", "qual"), "cigar"]}
    reader = AlignmentReader(path, reference=reference, cache=cache)
    # Each selection decodes the data series it needs, once.
    for fields, want in expected.items():
        assert snapshot(reader, fields) == want
    loads = cache.stats()["loads"]
    assert loads > 0
    for fields, want in expected.items():
        assert snapshot(reader, fields) == want
    assert cache.stats()["loads"] == loads
    clone = reader.clone()
    reader.close()
    try:
        assert snapshot(clone) == expected["all"]
        assert cache.stats()["loads"] == loads
        cache.clear()
        assert snapshot(clone) == expected["all"]
        assert cache.stats()["loads"] > loads
        clone.set_cache(None)
        loads = cache.stats()["loads"]
        assert snapshot(clone) == expected["all"]
        assert cache.stats()["loads"] == loads
    finally:
        clone.close()


@pytest.mark.parametrize("budget", [0, 1])
def test_cram_tiny_cache_and_failed_query_do_not_reuse_old_output(source, budget):
    path, reference = source
    cache = ow.cache.BlockCache(max_bytes=budget)
    with (
        AlignmentReader(path, reference=reference) as plain,
        AlignmentReader(path, reference=reference, cache=cache) as reader,
    ):
        for fields in [(), "all", ("seq", "qual"), "cigar", "all"]:
            assert snapshot(reader, fields) == snapshot(plain, fields)
        assert cache.stats()["resident_entries"] == 0
        with pytest.raises(ow.errors.UnknownChromError):
            reader.fetch(("missing", 0, 1), "all")
        assert snapshot(reader) == snapshot(plain)
        with pytest.raises(ValueError, match="field"):
            reader.fetch(("chr1", 0, 100), "unknown")
        empty = reader.fetch(("chr1", 20, 20), "all")
        assert len(empty) == 0
        assert empty.seqs.data.size == 0


def _read_child(reader, pipe):
    try:
        result = snapshot(reader)
        assert snapshot(reader) == result
        pipe.send(("ok", result))
    except BaseException as error:
        pipe.send(("error", type(error).__name__))
    finally:
        pipe.close()


@pytest.mark.parametrize("method", ["fork", "spawn"])
def test_cram_fork_rejects_and_spawn_reopens_without_cache(source, method):
    if method not in multiprocessing.get_all_start_methods():
        pytest.skip("start method unavailable")
    path, reference = source
    cache = ow.cache.BlockCache(max_bytes=16 << 20)
    with AlignmentReader(path, reference=reference, cache=cache) as reader:
        want = snapshot(reader)
        context = multiprocessing.get_context(method)
        parent, child = context.Pipe(duplex=False)
        process = context.Process(target=_read_child, args=(reader, child))
        process.start()
        child.close()
        try:
            assert parent.poll(30)
            got = parent.recv()
            assert got == (("error", "RuntimeError") if method == "fork" else ("ok", want))
        finally:
            process.join(5)
            if process.is_alive():
                process.terminate()
                process.join(5)
            parent.close()
        assert process.exitcode == 0
        reopened = pickle.loads(pickle.dumps(reader))
        try:
            assert snapshot(reopened) == want
        finally:
            reopened.close()


def test_cram_multireference_fields_and_cached_outputs(tmp_path):
    pysam = pytest.importorskip("pysam")
    # A reference-less file with both references in one slice exercises
    # per-reference decoded cache variants independently of a FASTA handle.
    source = tmp_path / "multi.cram"
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "a", "LN": 100}, {"SN": "b", "LN": 100}]}
    with pysam.AlignmentFile(
        str(source),
        "wc",
        header=header,
        format_options=[b"no_ref=1", b"multi_seq_per_slice=1"],
    ) as writer:
        for tid in range(2):
            record = pysam.AlignedSegment(writer.header)
            record.query_name = f"read{tid}"
            record.reference_id = tid
            record.reference_start = 10
            record.cigarstring = "4M"
            record.query_sequence = "ACGT"
            record.query_qualities = [20, 21, 22, 23]
            writer.write(record)
    ow.formats.Cram(source).build_index()
    cache = ow.cache.BlockCache(max_bytes=16 << 20)
    with AlignmentReader(source, cache=cache) as reader:
        for chrom in ["a", "b", "a", "b"]:
            region = (chrom, 0, 100)
            got = reader.fetch(region, "all")
            assert got.names.to_list() == [f"read{0 if chrom == 'a' else 1}"]
            assert got.seqs.to_list() == ["ACGT"]
            np.testing.assert_array_equal(got.quals.data, [20, 21, 22, 23])
        assert cache.stats()["hits"] > 0


def _itf8(data, at):
    byte = data[at]
    n = 1 if byte < 128 else 2 if byte < 192 else 3 if byte < 224 else 4 if byte < 240 else 5
    if n == 5:
        value = ((byte & 15) << 28) | (int.from_bytes(data[at + 1 : at + 4], "big") << 4)
        value |= data[at + 4] & 15
    else:
        value = byte & ((1 << (8 - n)) - 1)
        for b in data[at + 1 : at + n]:
            value = value << 8 | b
    return value, at + n


def _first_compression_header_crc(path):
    index = ow.formats.Crai.beside(ow.formats.Cram(path)).path.read_bytes()
    if index.startswith(b"\x1f\x8b"):
        index = gzip.decompress(index)
    at = int(index.splitlines()[0].split()[3]) + 4
    data = path.read_bytes()
    for _ in range(4):
        _, at = _itf8(data, at)
    for _ in range(2):
        byte = data[at]
        count = 0
        while count < 8 and byte & (128 >> count):
            count += 1
        at += count + 1
    _, at = _itf8(data, at)
    count, at = _itf8(data, at)
    for _ in range(count):
        _, at = _itf8(data, at)
    at += 4  # container header CRC
    assert data[at + 1] == 1  # compression-header block
    at += 2
    _, at = _itf8(data, at)
    compressed, at = _itf8(data, at)
    _, at = _itf8(data, at)
    return at + compressed


def test_cram_crc_failure_is_not_published_and_retry_recovers(source):
    path, reference = source
    original = path.read_bytes()
    corrupt = bytearray(original)
    corrupt[_first_compression_header_crc(path)] ^= 1
    cache = ow.cache.BlockCache(max_bytes=16 << 20)
    with AlignmentReader(path, reference=reference, cache=cache) as reader:
        expected = snapshot(reader)
        cache.clear()
        path.write_bytes(corrupt)
        try:
            with pytest.raises(ow.errors.InflateError, match="CRC32"):
                snapshot(reader)
        finally:
            path.write_bytes(original)
        assert snapshot(reader) == expected

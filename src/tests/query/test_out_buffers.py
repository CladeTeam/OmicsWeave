"""Caller storage is validated before writes and remains reusable after errors."""

import numpy as np
import pytest

from omicsweave.formats import BigWigWriter, FastaWriter, HicHeader, HicWriter
from omicsweave.query import ContactReader, SequenceReader, TrackReader
from omicsweave.types import ContactWindow


@pytest.fixture(params=["sequence", "signal", "bins", "contact"])
def dense_query(tmp_path, request):
    kind = request.param
    region = ("c", 0, 8)
    if kind == "sequence":
        path = tmp_path / "reference.fa"
        with FastaWriter(path) as writer:
            writer.add("c", b"ACGTACGT")
        reader = SequenceReader(path)
        expected = np.frombuffer(b"ACGTACGT", np.uint8)
        args = (region,)
        method, column = reader.fetch, "bases"
    elif kind in ("signal", "bins"):
        path = tmp_path / "signal.bw"
        with BigWigWriter(path, {"c": 8}) as writer:
            writer.add("c", [0, 4], [4, 8], [2.0, 3.0])
        reader = TrackReader(path)
        if kind == "signal":
            method, args = reader.signal, (region,)
            expected = np.array([2.0] * 4 + [3.0] * 4, np.float32)
        else:
            method, args = reader.bins, (region, 4, "sum")
            expected = np.array([[8.0, 12.0]], np.float32)
        column = "values"
    else:
        path = tmp_path / "contacts.hic"
        with HicWriter(path, HicHeader("buffers", {"a": 8, "b": 8}, (4,))) as writer:
            writer.write_matrix("a", "b", 4, [0, 1], [0, 1], [2.0, 3.0])
        reader = ContactReader(path)
        method = reader.matrix
        args = (ContactWindow(("a", 0, 8), ("b", 0, 8), 4),)
        expected = np.array([[2.0, 0.0], [0.0, 3.0]], np.float32)
        column = "values"
    try:
        yield reader, method, args, column, expected
    finally:
        reader.close()


@pytest.mark.parametrize("invalid", ["dtype", "strides", "readonly", "shape"])
def test_invalid_out_preserves_storage_and_reader_can_retry(dense_query, invalid):
    _, query, args, column, expected = dense_query
    storage = np.full(expected.shape, 17, expected.dtype)
    if invalid == "dtype":
        storage = storage.astype(np.float64)
        output = storage
    elif invalid == "strides":
        storage = np.full((*expected.shape[:-1], expected.shape[-1] * 2), 17, expected.dtype)
        output = storage[..., ::2]
        assert not output.flags.c_contiguous
    elif invalid == "readonly":
        storage.flags.writeable = False
        output = storage
    else:
        output = storage.reshape((expected.size, 1) if expected.ndim == 1 else (expected.size,))
    before = storage.copy()
    with pytest.raises((TypeError, ValueError)):
        query(*args, out=output)
    np.testing.assert_array_equal(storage, before)

    reusable = np.full(expected.shape, 17, expected.dtype)
    assert getattr(query(*args, out=reusable), column) is reusable
    np.testing.assert_array_equal(reusable, expected)


def test_closed_reader_rejects_out_without_mutating_it(dense_query):
    reader, query, args, _, expected = dense_query
    reader.close()
    reader.close()
    output = np.full(expected.shape, 17, expected.dtype)
    with pytest.raises(ValueError, match="closed"):
        query(*args, out=output)
    np.testing.assert_array_equal(output, 17)

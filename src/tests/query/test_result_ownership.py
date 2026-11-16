import numpy as np
import pytest

from omicsweave.cache import BlockCache
from omicsweave.formats import BigWigWriter, FastaWriter, HicHeader, HicWriter
from omicsweave.query import ContactReader, SequenceReader, TrackReader
from omicsweave.types import ContactWindow


@pytest.mark.parametrize("cached", [False, True])
@pytest.mark.parametrize("length", [0, 8, 8192])
def test_sequence_results_survive_queries_cache_clear_and_close(tmp_path, cached, length):
    path = tmp_path / "reference.fa.gz"
    bases = b"ACGT" * 4096
    with FastaWriter(path) as writer:
        writer.add("c", bases)
    cache = BlockCache(max_bytes=1 << 20) if cached else None
    with SequenceReader(path, cache=cache) as reader:
        first = reader.fetch(("c", 0, length))
        second = reader.fetch(("c", 0, length))
        second.bases[:] = ord("N")
        reader.fetch(("c", 1, 9))
        if cache is not None:
            cache.clear()
    assert first.bases.dtype == np.uint8
    assert first.bases.shape == (length,)
    assert first.bases.tobytes() == bases[:length]
    assert second.bases.tobytes() == b"N" * length


@pytest.mark.parametrize("cached", [False, True])
def test_signal_and_bins_own_results_while_out_shares_storage(tmp_path, cached):
    path = tmp_path / "signal.bw"
    with BigWigWriter(path, {"c": 16}) as writer:
        writer.add("c", [0, 8], [4, 12], [2.0, 3.0])
    cache = BlockCache(max_bytes=1 << 20) if cached else None
    with TrackReader(path, cache=cache) as reader:
        signal = reader.signal(("c", 0, 16))
        bins = reader.bins(("c", 0, 16), 4, "sum")
        reader.signal(("c", 0, 16)).values[:] = -1
        reader.bins(("c", 0, 16), 4, "sum").values[:] = -1
        storage = np.empty(16, np.float32)
        shared = reader.signal(("c", 0, 16), out=storage)
        assert shared.values is storage
        storage[:] = -2
        if cache is not None:
            cache.clear()
    expected = [2.0] * 4 + [np.nan] * 4 + [3.0] * 4 + [np.nan] * 4
    assert signal.values.dtype == bins.values.dtype == np.float32
    np.testing.assert_array_equal(signal.values, expected)
    np.testing.assert_array_equal(bins.values, [[8.0, np.nan, 12.0, np.nan]])
    np.testing.assert_array_equal(shared.values, -2)


@pytest.mark.parametrize("cached", [False, True])
def test_contact_results_survive_queries_cache_clear_and_close(tmp_path, cached):
    path = tmp_path / "contacts.hic"
    with HicWriter(path, HicHeader("test", {"a": 20, "b": 20}, (10,))) as writer:
        writer.write_matrix("a", "b", 10, [0, 1], [1, 0], [2.0, 3.0])
    cache = BlockCache(max_bytes=1 << 20) if cached else None
    window = ContactWindow(("a", 0, 20), ("b", 0, 20), 10)
    with ContactReader(path, cache=cache) as reader:
        first = reader.matrix(window)
        reader.matrix(window).values[:] = -1
        reader.pixels(window)
        if cache is not None:
            cache.clear()
    assert first.values.dtype == np.float32
    np.testing.assert_array_equal(first.values, [[0.0, 2.0], [3.0, 0.0]])

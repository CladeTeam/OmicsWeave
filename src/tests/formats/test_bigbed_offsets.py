import numpy as np
import pytest

from omicsweave import _native
from omicsweave.query import IntervalReader


def native_writer(path):
    return _native.BigBedWriter(str(path), [("c", 100)], 4, 4, None, 2, 3, False, 6, 0, None)


@pytest.mark.parametrize(
    "offsets",
    [
        [-1, 2, 3],
        [0, -1, 3],
        [0, 2, 1],
        [0, 2, 4],
        [0, 2, np.iinfo(np.int64).max],
    ],
)
def test_invalid_offsets_validate_whole_batch_before_writing(tmp_path, offsets):
    path = tmp_path / "offsets.bb"
    writer = native_writer(path)
    starts = np.array([1, 3], np.uint32)
    ends = starts + 1
    data = np.frombuffer(b"abc", np.uint8)
    with pytest.raises(ValueError, match="annotation offsets are invalid"):
        writer.add(0, starts, ends, np.array(offsets, np.int64), data)
    writer.add(0, starts, ends, np.array([0, 2, 3], np.int64), data)
    writer.close(False)
    with IntervalReader(path) as reader:
        result = reader.fetch(("c", 0, 100))
        np.testing.assert_array_equal(result.start, starts)
        assert result.annotation.to_list() == ["ab", "c"]


def test_offset_subranges_empty_batches_and_input_lifetime(tmp_path):
    path = tmp_path / "subranges.bb"
    writer = native_writer(path)
    empty = np.array([], np.uint32)
    writer.add(0, empty, empty, np.array([99], np.int64), np.array([], np.uint8))
    starts = np.array([1, 3, 5], np.uint32)
    ends = starts + 1
    offsets = np.array([1, 1, 3, 4], np.int64)
    data = np.frombuffer(b"xabcy", np.uint8).copy()
    writer.add(0, starts, ends, offsets, data)
    starts[:] = 99
    ends[:] = 99
    offsets[:] = 0
    data[:] = ord("z")
    del starts, ends, offsets, data
    writer.close(False)
    assert writer.close(False) is None
    with IntervalReader(path) as reader:
        result = reader.fetch(("c", 0, 100))
        np.testing.assert_array_equal(result.start, [1, 3, 5])
        np.testing.assert_array_equal(result.end, [2, 4, 6])
        assert result.annotation.to_list() == ["", "ab", "c"]
    with pytest.raises(ValueError, match="closed file"):
        writer.add(0, empty, empty, np.array([0], np.int64), np.array([], np.uint8))


def test_offset_column_layout_and_length_errors_leave_writer_usable(tmp_path):
    path = tmp_path / "layout.bb"
    writer = native_writer(path)
    starts = np.array([1, 3], np.uint32)
    ends = starts + 1
    data = np.frombuffer(b"abc", np.uint8)
    with pytest.raises(ValueError, match="differ in length"):
        writer.add(0, starts, ends, np.array([0, 3], np.int64), data)
    offsets = np.array([0, 99, 2, 99, 3], np.int64)[::2]
    with pytest.raises(ValueError, match="contiguous"):
        writer.add(0, starts, ends, offsets, data)
    writer.add(0, starts, ends, np.array([0, 2, 3], np.int64), data)
    writer.close(False)
    with IntervalReader(path) as reader:
        assert reader.fetch(("c", 0, 100)).annotation.to_list() == ["ab", "c"]

import dataclasses

import numpy as np
import pytest

from omicsweave import errors
from omicsweave.types import BinUnit, ContactWindow, Region


@pytest.mark.parametrize("resolution", [True, False, 1.5, float("nan"), "10", None])
def test_contact_window_rejects_noninteger_resolution(resolution):
    with pytest.raises(TypeError, match="resolution must be an integer"):
        ContactWindow(("a", 0, 20), ("b", 0, 30), resolution)


@pytest.mark.parametrize("resolution", [0, -1, 2**32])
def test_contact_window_rejects_resolution_outside_u32(resolution):
    with pytest.raises(errors.InvalidError):
        ContactWindow(("a", 0, 20), ("b", 0, 30), resolution)


@pytest.mark.parametrize("resolution", [10, np.int32(10), np.uint64(10)])
def test_contact_window_normalizes_integer_resolution(resolution):
    window = ContactWindow(("a", 0, 20), ("b", 0, 30), resolution)
    assert type(window.resolution) is int
    assert window.shape == (2, 3)
    assert window.row_bins == range(2)
    assert window.col_bins == range(3)


def test_region_and_window_keys():
    r = Region("chr1", np.int64(5), 9)
    assert (r.start, len(r), str(r)) == (5, 4, "chr1:5-9")
    assert dataclasses.replace(r, end=20).end == 20
    assert hash(r) == hash(Region("chr1", 5, 9))
    with pytest.raises(ValueError):
        Region("chr1", -1, 3)
    with pytest.raises(TypeError):
        Region("chr1", 1.5, 3)
    w = ContactWindow(("a", 0, 2500), ("b", 999, 1001), 1000, "FRAG")
    assert w.unit is BinUnit.FRAG
    assert w.rows == Region("a", 0, 2500)
    assert w.shape == (3, 2)
    assert ContactWindow(("a", 7, 7), ("a", 0, 1), 5).shape == (0, 1)

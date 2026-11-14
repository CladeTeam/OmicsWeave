import numpy as np
import pytest

from omicsweave import errors
from omicsweave.types import ContactWindow


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

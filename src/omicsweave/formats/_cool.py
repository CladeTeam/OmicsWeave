"""cooler (HDF5): one contact map (``.cool``), or the same map at several
resolutions (``.mcool``)."""

from __future__ import annotations

from ._base import File

__all__ = ["Cool", "Mcool"]


class Cool(File):
    """A single-resolution cooler file (``.cool``)."""

    __slots__ = ()
    format = "cool"


class Mcool(File):
    """A multi-resolution cooler file (``.mcool``): one cooler per
    resolution under ``/resolutions``."""

    __slots__ = ()
    format = "mcool"

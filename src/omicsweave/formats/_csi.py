"""Coordinate-sorted index (``.csi``): where the records of each genomic bin
start in a BGZF file (BCF, BAM)."""

from __future__ import annotations

from ._base import File, PathLike, companion_path

__all__ = ["Csi"]


class Csi(File):
    """A ``.csi`` index of a sorted BCF or BAM file, as htslib writes it."""

    __slots__ = ()
    format = "csi"
    suffix = ".csi"

    @classmethod
    def beside(cls, indexed: PathLike) -> Csi:
        """The index named after ``indexed``: ``<indexed>.csi``."""
        return cls(companion_path(indexed, cls.suffix))

"""bgzip index (``.gzi``): compressed and uncompressed offsets of every BGZF block."""

from __future__ import annotations

from ._base import File, PathLike, companion_path

__all__ = ["Gzi"]


class Gzi(File):
    """A ``.gzi`` index of a BGZF file, as ``bgzip -i`` writes it.

    Bgzip FASTA needs one, besides its ``.fai``, to be read by region.
    """

    __slots__ = ()
    format = "gzi"
    suffix = ".gzi"

    @classmethod
    def beside(cls, bgzf: PathLike) -> Gzi:
        """The index named after ``bgzf``: ``<bgzf>.gzi``."""
        return cls(companion_path(bgzf, cls.suffix))

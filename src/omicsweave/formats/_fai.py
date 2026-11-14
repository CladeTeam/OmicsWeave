"""FASTA index (``.fai``): where each sequence's lines start, as samtools writes it."""

from __future__ import annotations

from ._base import File, PathLike, companion_path

__all__ = ["Fai"]


class Fai(File):
    """A ``.fai`` FASTA index: name, length, offset and line layout of each sequence.

    Plain and bgzip FASTA both need one to be read by region.
    """

    __slots__ = ()
    format = "fai"
    suffix = ".fai"

    @classmethod
    def beside(cls, fasta: PathLike) -> Fai:
        """The index named after ``fasta``: ``<fasta>.fai``."""
        return cls(companion_path(fasta, cls.suffix))

"""BAM index (``.bai``): where the records of each genomic bin start in a
BAM file."""

from __future__ import annotations

from ._base import File, PathLike, companion_path

__all__ = ["Bai"]


class Bai(File):
    """A ``.bai`` index of a sorted BAM file, as ``samtools index`` writes it."""

    __slots__ = ()
    format = "bai"
    suffix = ".bai"

    @classmethod
    def beside(cls, bam: PathLike) -> Bai:
        """The index named after ``bam``: ``<bam>.bai``."""
        return cls(companion_path(bam, cls.suffix))

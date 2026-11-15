"""CRAM index (``.crai``): where each slice of a CRAM file lies, as
gzip-compressed text."""

from __future__ import annotations

from ._base import File, PathLike, companion_path

__all__ = ["Crai"]


class Crai(File):
    """A ``.crai`` index of a sorted CRAM file, as ``samtools index`` writes
    it."""

    __slots__ = ()
    format = "crai"
    suffix = ".crai"

    @classmethod
    def beside(cls, cram: PathLike) -> Crai:
        """The index named after ``cram``: ``<cram>.crai``."""
        return cls(companion_path(cram, cls.suffix))

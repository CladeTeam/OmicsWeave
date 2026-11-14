"""UCSC 2bit (``.2bit``): packed bases with N and soft-mask blocks."""

from __future__ import annotations

from ._base import File

__all__ = ["TwoBit"]


class TwoBit(File):
    """A ``.2bit`` file; it indexes itself, so it is read by region as it is."""

    __slots__ = ()
    format = "2bit"

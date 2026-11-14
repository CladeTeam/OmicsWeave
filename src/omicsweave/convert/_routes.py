"""Which conversions there are."""

from __future__ import annotations

from .. import _native
from ..formats import FORMATS, File


def routes() -> tuple[tuple[type[File], type[File]], ...]:
    """Every (source, target) pair of formats there is a conversion for."""
    return tuple((FORMATS[a], FORMATS[b]) for a, b in _native.conversions())

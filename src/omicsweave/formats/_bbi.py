"""The BBI container that bigWig and bigBed share: how its blocks, index and
zoom levels are written."""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["BbiOptions"]


@dataclass(frozen=True, slots=True)
class BbiOptions:
    """How a bigWig or bigBed file is written (UCSC defaults)."""

    items_per_slot: int | None = None
    """Items per data block; ``None``: 1024 for bigWig, 512 for bigBed."""
    block_size: int = 256
    """Children per index node."""
    compress: bool = True
    level: int = 6
    """Compression level, 1 to 12."""
    zoom_levels: int = 10
    """Most zoom levels derived."""
    zoom_reductions: tuple[int, ...] | None = None
    """Explicit zoom reductions instead of derived ones."""


def native_options(options: BbiOptions) -> tuple[object, ...]:
    reductions = None if options.zoom_reductions is None else list(options.zoom_reductions)
    return (
        options.items_per_slot,
        options.block_size,
        options.compress,
        options.level,
        options.zoom_levels,
        reductions,
    )

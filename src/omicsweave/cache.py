"""A store of decoded blocks that readers share.

Give readers one :class:`BlockCache` (``cache=``) and the blocks one query
decodes serve the next, across readers and their clones: a sliding window,
a sequence and its signal at the same loci, many workers' threads. Which
blocks a query needs, and when they are loaded or dropped, is the native
layer's business.

The budget covers admitted payload and prepared indexes, including blocks a
running query holds; it excludes decoder scratch, in-flight allocations and
allocator overhead, so it is not a limit on process memory.
"""

from __future__ import annotations

import operator

from . import _native

__all__ = ["BlockCache"]


def _size(value: int, name: str) -> int:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an integer")
    number = operator.index(value)
    if number < 0:
        raise ValueError(f"{name} must be nonnegative")
    return number


class BlockCache:
    """Decoded blocks shared by the readers given it, up to ``max_bytes``
    and ``max_entries`` (FIFO admission; 0 bytes caches nothing).

    The files a cache serves must stay unchanged while readers use it. A
    cache is process-local: it cannot be pickled or used after fork, so
    create one in each worker process.
    """

    __slots__ = ("_native",)

    def __init__(self, *, max_bytes: int, max_entries: int = 8192) -> None:
        self._native = _native.BlockCache(
            _size(max_bytes, "max_bytes"), _size(max_entries, "max_entries")
        )

    def stats(self) -> dict[str, int]:
        """Hits, misses, loads, evictions, what is resident, and the budgets."""
        return self._native.stats()

    def clear(self) -> None:
        """Drops every entry; blocks a running query holds stay valid."""
        self._native.clear()

    def resize(self, *, max_bytes: int, max_entries: int) -> None:
        """Sets new budgets; blocks a running query holds stay valid."""
        self._native.resize(_size(max_bytes, "max_bytes"), _size(max_entries, "max_entries"))

    def __reduce__(self) -> object:
        raise TypeError("BlockCache is process-local; create one in each worker")


def native(cache: BlockCache | None) -> _native.BlockCache | None:
    if cache is not None and not isinstance(cache, BlockCache):
        raise TypeError("cache must be a BlockCache or None")
    return None if cache is None else cache._native

"""Explicit, process-local ownership of reusable decoded blocks.

Cache budgets cover admitted payload capacity and prepared indexes, including
live borrowed blocks. They exclude decoder scratch, in-flight allocations and
allocator/map overhead, and therefore are not process RSS limits.
"""

from __future__ import annotations

import operator

from . import _native

BlockKey = _native.BlockKey
CacheNamespace = _native.CacheNamespace
__all__ = ["BlockCache", "BlockKey", "CacheNamespace"]


def _size(value: int, name: str) -> int:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an integer")
    number = operator.index(value)
    if number < 0:
        raise ValueError(f"{name} must be nonnegative")
    return number


class BlockCache:
    """Caller-owned block storage, shared explicitly among readers.

    Inputs and indexes must remain unchanged during each reader session.
    Independent opens have separate namespaces. Cloned readers share theirs.
    FIFO is the convenience admission policy; explicit eviction lets a loader
    choose residency. Objects and their contents cannot be pickled or reused
    after fork. Create a cache separately in each worker process.
    """

    __slots__ = ("_native",)

    def __init__(self, *, max_bytes: int, max_entries: int = 8192) -> None:
        self._native = _native.BlockCache(
            _size(max_bytes, "max_bytes"), _size(max_entries, "max_entries")
        )

    def stats(self) -> dict[str, int]:
        """Residency, live leases, loads and current admission budgets."""
        return self._native.stats()

    def keys(self) -> tuple[BlockKey, ...]:
        """The entries currently available for lookup."""
        return tuple(self._native.keys())

    def evict(self, key: BlockKey) -> bool:
        """Remove an entry without invalidating an active query's lease."""
        return self._native.evict(key)

    def clear(self, *, source: CacheNamespace | None = None) -> None:
        """Clear all entries, or one source; in-flight loads cannot republish."""
        self._native.clear(source)

    def resize(self, *, max_bytes: int, max_entries: int) -> None:
        """Set new admission budgets; active leases remain valid and charged."""
        self._native.resize(_size(max_bytes, "max_bytes"), _size(max_entries, "max_entries"))

    def __reduce__(self) -> object:
        raise TypeError("BlockCache is process-local; create one in each worker")


def native(cache: BlockCache | None) -> _native.BlockCache | None:
    if cache is not None and not isinstance(cache, BlockCache):
        raise TypeError("cache must be a BlockCache or None")
    return None if cache is None else cache._native

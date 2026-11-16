"""Defaults of the native layer, read once."""

from __future__ import annotations

from collections.abc import Callable
from functools import cache
from typing import Any

from .. import _native

__all__ = ["default"]


@cache
def _all() -> dict[str, dict[str, Any]]:
    return _native.motif_defaults()


def default(group: str, name: str) -> Callable[[], Any]:
    """A dataclass ``default_factory`` for the native default ``group.name``."""
    return lambda: _all()[group][name]

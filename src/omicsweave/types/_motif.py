"""Motif models as FFM stores them: an alphabet, and per motif a window of
positions scored by local factors."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Literal

import numpy as np

from .. import _native

__all__ = ["Alphabet", "AnnotationValue", "Factor", "FactorKind", "Gaps", "Motif"]

FactorKind = Literal["joint", "conditional", "potential"]
"""How a factor's table reads: ``"joint"`` frequencies of the symbol
combination (normalized over the table), ``"conditional"`` frequencies of
the last position's symbol given the earlier ones (normalized per context),
``"potential"`` unnormalized weights (score ``ln(value)``; 0 forbids)."""

AnnotationValue = str | int | float
"""An annotation value, typed by the collection's schema."""


@dataclass(frozen=True, slots=True, eq=False)
class Alphabet:
    """An ordered set of one-character symbols; symbol ``i`` has code ``i``.

    ``complement`` names the symbol paired with each symbol (``"TGCA"`` for
    ``"ACGT"``); without it only the forward strand can be scanned.
    ``aliases`` maps other spellings to symbols (``{"a": "A"}``). Any other
    byte is outside the alphabet.
    """

    name: str
    symbols: str
    complement: str | None = None
    aliases: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "aliases", MappingProxyType(dict(self.aliases)))

    @classmethod
    def dna(cls) -> Alphabet:
        """``ACGT``, lowercase accepted, with the Watson-Crick complement."""
        return cls._preset("dna")

    @classmethod
    def rna(cls) -> Alphabet:
        """``ACGU``, lowercase accepted, with the Watson-Crick complement."""
        return cls._preset("rna")

    @classmethod
    def protein(cls) -> Alphabet:
        """The twenty standard amino acids, lowercase accepted."""
        return cls._preset("protein")

    @classmethod
    def _preset(cls, name: str) -> Alphabet:
        return cls(*_native.alphabet_preset(name))

    def __len__(self) -> int:
        return len(self.symbols)

    @property
    def lut(self) -> np.ndarray:
        """Code of every byte (``uint8``, 256 entries); 255 outside the alphabet."""
        lut = np.full(256, 255, dtype=np.uint8)
        for code, s in enumerate(self.symbols):
            lut[ord(s)] = code
        for byte, s in self.aliases.items():
            lut[ord(byte)] = self.symbols.index(s)
        return lut

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Alphabet):
            return NotImplemented
        return (self.name, self.symbols, self.complement, dict(self.aliases)) == (
            other.name,
            other.symbols,
            other.complement,
            dict(other.aliases),
        )

    def __hash__(self) -> int:
        return hash((self.name, self.symbols, self.complement))

    def __repr__(self) -> str:
        comp = f", complement={self.complement!r}" if self.complement else ""
        return f"Alphabet({self.name!r}, {self.symbols!r}{comp})"


@dataclass(frozen=True, slots=True, eq=False)
class Factor:
    """A table over the local codes at a few positions of the motif window.

    ``positions`` strictly increase; ``table`` (``float64``) has
    ``arity ** len(positions)`` cells, the last position varying fastest,
    where ``arity`` is the motif's (its support, plus one rest category when
    the support is a strict subset of the alphabet).
    """

    kind: FactorKind
    positions: tuple[int, ...]
    table: np.ndarray

    def __post_init__(self) -> None:
        if self.kind not in ("joint", "conditional", "potential"):
            raise ValueError(
                f"factor kind must be joint, conditional or potential, not {self.kind!r}"
            )
        object.__setattr__(self, "positions", tuple(int(p) for p in self.positions))
        object.__setattr__(self, "table", np.ascontiguousarray(self.table, dtype=np.float64))

    @property
    def order(self) -> int:
        return len(self.positions)


@dataclass(frozen=True, slots=True, eq=False)
class Gaps:
    """Gap counts of a gapped profile (an alignment of ``sites`` sequences):
    ``deletions`` per column (sequences missing it, ``width`` values) and
    ``insertions`` per junction between consecutive columns (residues
    inserted there over all sequences, ``width - 1`` values)."""

    deletions: np.ndarray
    insertions: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "deletions", np.ascontiguousarray(self.deletions, dtype=np.float64)
        )
        object.__setattr__(
            self, "insertions", np.ascontiguousarray(self.insertions, dtype=np.float64)
        )


@dataclass(frozen=True, slots=True, eq=False)
class Motif:
    """One motif: its ID, its model and its annotation.

    The model scores a window of ``width`` positions with ``factors``.
    ``support`` lists the codes the motif models individually (ascending);
    every other symbol falls into one rest category. ``None`` supports the
    whole alphabet. ``sites`` counts the observations behind the
    frequencies, when known. ``gaps`` makes a gapped profile, whose factors
    are one single-position joint factor per column. The annotation holds
    the values of the collection's schema keys; missing ones are absent.
    """

    id: str
    width: int
    factors: tuple[Factor, ...]
    support: bytes | None = None
    sites: float | None = None
    gaps: Gaps | None = None
    annotation: Mapping[str, AnnotationValue] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "factors", tuple(self.factors))
        if self.support is not None:
            object.__setattr__(self, "support", bytes(self.support))
        object.__setattr__(self, "annotation", MappingProxyType(dict(self.annotation)))

    @classmethod
    def from_pfm(
        cls,
        id: str,
        rows: np.ndarray | Iterable[Iterable[float]],
        *,
        sites: float | None = None,
        annotation: Mapping[str, AnnotationValue] | None = None,
    ) -> Motif:
        """A position frequency matrix: one row of counts or frequencies per
        position, one column per symbol of the alphabet."""
        matrix = np.asarray(rows, dtype=np.float64)
        if matrix.ndim != 2:
            raise ValueError("a frequency matrix has one row per position")
        factors = tuple(Factor("joint", (i,), row) for i, row in enumerate(matrix))
        return cls(id, len(matrix), factors, sites=sites, annotation=annotation or {})

    @property
    def matrix(self) -> np.ndarray | None:
        """The ``(width, alphabet)`` frequency matrix when the motif is one
        single-position joint factor per position over the whole alphabet,
        in position order (a gapped profile's residue counts); ``None``
        otherwise."""
        if self.support is not None or len(self.factors) != self.width:
            return None
        for i, f in enumerate(self.factors):
            if f.kind != "joint" or f.positions != (i,):
                return None
        if not self.factors:
            return None
        return np.stack([f.table for f in self.factors])

    @property
    def gapped(self) -> bool:
        return self.gaps is not None

    def __repr__(self) -> str:
        kind = "gapped profile" if self.gaps is not None else f"{len(self.factors)} factors"
        return f"Motif({self.id!r}, width={self.width}, {kind})"

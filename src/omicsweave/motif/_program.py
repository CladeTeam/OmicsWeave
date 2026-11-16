"""Compiled motif programs: every motif's scores, p-value tables and scan
tables, ready to scan sequences once for all of them."""

from __future__ import annotations

import os
from collections.abc import Generator, Iterable
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any, overload

import numpy as np

from .. import _native
from ..formats import Ffm
from ..query import MotifReader
from ..types import Alphabet, Motif
from ._background import Background
from ._compiled import CompiledMotif, NullDistribution, Term
from ._reduce import Reduction, build_result
from ._scoring import Scoring, Strands
from ._sequences import Sequences, as_strings

__all__ = ["MotifSource", "Part", "Program", "combine", "compile"]

MotifSource = MotifReader | Ffm | str | os.PathLike[str] | Iterable[Motif]
"""Motifs to compile: an open :class:`~omicsweave.query.MotifReader` (its
alphabet is used and no Python objects are built), an FFM file, or
:class:`~omicsweave.types.Motif` values with an explicit alphabet."""


@contextmanager
def opened(
    motifs: MotifSource, alphabet: Alphabet | None
) -> Generator[tuple[MotifReader | list[Motif], Alphabet]]:
    """``motifs`` ready to compile, with their alphabet; a file stays open
    for the ``with`` block."""
    if isinstance(motifs, str | os.PathLike):
        with MotifReader(motifs) as reader, opened(reader, alphabet) as out:
            yield out
        return
    if isinstance(motifs, MotifReader):
        if alphabet is not None and alphabet != motifs.alphabet:
            raise ValueError("alphabet differs from the reader's")
        yield motifs, motifs.alphabet
        return
    if alphabet is None:
        raise TypeError("motifs given as values need an alphabet")
    values = list(motifs)
    for m in values:
        if not isinstance(m, Motif):
            raise TypeError(f"expected Motif values, not {type(m).__name__}")
    yield values, alphabet


def native_source(motifs: MotifReader | list[Motif]) -> Any:
    """What the native compilers take for ``motifs``."""
    return motifs._native if isinstance(motifs, MotifReader) else motifs


@dataclass(frozen=True, slots=True, eq=False)
class Part:
    """Motifs compiled under one policy: motifs ``motifs`` of the program."""

    background: Background
    scoring: Scoring
    strands: Strands
    motifs: range


class Program:
    """Compiled motifs. Scanning a program reads each sequence once for all
    its motifs; what is kept is decided by the reductions passed to
    :meth:`scan`. A program is immutable and can be shared across threads;
    scans release the GIL.
    """

    __slots__ = ("_native", "alphabet", "parts")

    alphabet: Alphabet
    parts: tuple[Part, ...]

    def __init__(self, native: Any, alphabet: Alphabet, parts: tuple[Part, ...]) -> None:
        self._native = native
        self.alphabet = alphabet
        self.parts = parts

    def __len__(self) -> int:
        return self._native.__len__()

    @property
    def widths(self) -> np.ndarray:
        """Each motif's width (``uint32``)."""
        return self._native.widths()

    def part_of(self, motif: int) -> int:
        """The part motif ``motif`` belongs to."""
        for i, p in enumerate(self.parts):
            if motif in p.motifs:
                return i
        raise IndexError("motif index out of range")

    def model(self, motif: int, strand: str = "+") -> CompiledMotif:
        """Motif ``motif`` as compiled for ``strand``: its terms, their
        integer form and its exact null distribution."""
        d = self._native.model(motif, strand)
        part = self.parts[self.part_of(motif)]
        terms = tuple(Term(tuple(p.tolist()), t) for p, t in d["terms"])
        int_terms = tuple(Term(tuple(p.tolist()), t) for p, t in d["int_terms"])
        null = NullDistribution(d["pmf"], d["survival"], d["forbidden"])
        return CompiledMotif(
            motif,
            strand,
            d["width"],
            terms,
            int_terms,
            d["scale"],
            d["offset"],
            d["max"],
            null,
            part.scoring.rounding,
        )

    @overload
    def scan(
        self,
        sequences: Sequences,
        reduction: Reduction,
        /,
        *,
        threads: int | None = ...,
        tile_len: int | None = ...,
    ) -> Any: ...
    @overload
    def scan(
        self,
        sequences: Sequences,
        reduction: Reduction,
        /,
        *reductions: Reduction,
        threads: int | None = ...,
        tile_len: int | None = ...,
    ) -> tuple[Any, ...]: ...
    def scan(
        self,
        sequences: Sequences,
        *reductions: Reduction,
        threads: int | None = None,
        tile_len: int | None = None,
    ) -> Any:
        """Scans ``sequences`` once and returns each reduction's result (one
        result for one reduction, else a tuple in order). Results do not
        depend on ``threads`` (``None``: the available parallelism) or
        ``tile_len`` (letters per work unit)."""
        if not reductions:
            raise TypeError("scan needs at least one reduction")
        for r in reductions:
            if not isinstance(r, Reduction):
                raise TypeError(f"expected a reduction, not {type(r).__name__}")
        s = as_strings(sequences)
        outputs = self._native.scan(s.offsets, s.data, list(reductions), threads, tile_len)
        for result in outputs:
            if result["kind"] == "scores":
                result["lengths"] = np.diff(s.offsets)
        results = tuple(build_result(r) for r in outputs)
        return results[0] if len(results) == 1 else results

    def __repr__(self) -> str:
        parts = f"{len(self.parts)} part{'s' * (len(self.parts) != 1)}"
        return f"Program({len(self)} motifs, {parts}, {self.alphabet.name})"


def compile(
    motifs: MotifSource,
    background: Background,
    scoring: Scoring | None = None,
    *,
    strands: Strands = "both",
    alphabet: Alphabet | None = None,
    threads: int | None = None,
) -> Program:
    """Compiles ``motifs`` against ``background`` under ``scoring``
    (:meth:`Scoring.fimo` by default) for scanning ``strands``: frequencies
    are smoothed and contrasted with the background, scores quantized,
    exact p-value tables computed and scan tables built. Motifs compile in
    parallel on ``threads`` (``None``: the available parallelism)."""
    scoring = Scoring.fimo() if scoring is None else scoring
    if strands not in ("both", "+", "-"):
        raise ValueError(f"strands must be both, + or -, not {strands!r}")
    with opened(motifs, alphabet) as (source, alphabet):
        native = _native.MotifProgram(
            native_source(source), alphabet, list(background.kmers), scoring, strands, threads
        )
    part = Part(background, scoring, strands, range(len(native)))
    return Program(native, alphabet, (part,))


def combine(*programs: Program) -> Program:
    """One program scanning every motif of ``programs`` (which share an
    alphabet) in one pass: part ``i`` of the result is the ``i``-th part met
    and motif numbers follow the same order. Reductions choose a part with
    ``part=``."""
    if not programs:
        raise TypeError("combine needs at least one program")
    alphabet = programs[0].alphabet
    if any(p.alphabet != alphabet for p in programs):
        raise ValueError("combined programs must share an alphabet")
    native = _native.MotifProgram.combine([p._native for p in programs])
    parts: list[Part] = []
    base = 0
    for p in programs:
        for part in p.parts:
            r = part.motifs
            parts.append(
                Part(
                    part.background,
                    part.scoring,
                    part.strands,
                    range(r.start + base, r.stop + base),
                )
            )
        base += len(p)
    return Program(native, alphabet, tuple(parts))

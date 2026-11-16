"""Gapped profiles (GLAM2 motifs) aligned locally to sequences."""

from __future__ import annotations

import dataclasses
from collections.abc import Sequence as Seq
from dataclasses import dataclass
from typing import Any

import numpy as np

from .. import _native
from ..types import Alphabet, StringVector
from ._background import Background
from ._program import MotifSource, native_source, opened
from ._scoring import Strands
from ._sequences import Sequences, as_strings

__all__ = [
    "DirichletMixture",
    "ProfileAlignments",
    "ProfileProgram",
    "ProfileScoring",
    "compile_profiles",
]


@dataclass(frozen=True, slots=True, eq=False)
class DirichletMixture:
    """A Dirichlet mixture prior over a column's residue probabilities:
    component ``c`` has weight ``weights[c]`` and pseudocounts
    ``alphas[c]`` (one per symbol)."""

    weights: np.ndarray
    alphas: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "weights", np.ascontiguousarray(self.weights, dtype=np.float64))
        object.__setattr__(
            self, "alphas", np.ascontiguousarray(np.atleast_2d(self.alphas), dtype=np.float64)
        )

    @classmethod
    def glam2(cls, alphabet: Alphabet) -> DirichletMixture:
        """GLAM2's built-in prior for the alphabet: its mixture for amino
        acids, one uniform component otherwise."""
        return cls(*_native.dirichlet_glam2(alphabet))

    @classmethod
    def uniform(cls, alphabet: Alphabet) -> DirichletMixture:
        """One component with every pseudocount 1."""
        return cls(*_native.dirichlet_uniform(len(alphabet)))

    @classmethod
    def parse_glam2(cls, text: str, alphabet: Alphabet) -> DirichletMixture:
        """A mixture file as ``glam2scan -d`` reads it."""
        return cls(*_native.dirichlet_parse_glam2(text, len(alphabet)))


@dataclass(frozen=True, slots=True, eq=False)
class ProfileScoring:
    """How gap and residue counts become alignment scores: beta
    pseudocounts for a column being deleted or present and for a residue
    inserted after a column or not (``glam2scan -D -E -I -J``), and the
    residue prior (``-d``)."""

    delete_pseudo: float
    no_delete_pseudo: float
    insert_pseudo: float
    no_insert_pseudo: float
    prior: DirichletMixture

    @classmethod
    def glam2(cls, alphabet: Alphabet) -> ProfileScoring:
        """``glam2scan``'s defaults for the alphabet."""
        d, e, i, j, (weights, alphas) = _native.profile_scoring_glam2(alphabet)
        return cls(d, e, i, j, DirichletMixture(weights, alphas))

    def replace(self, **changes: Any) -> ProfileScoring:
        """A copy with ``changes`` applied."""
        return dataclasses.replace(self, **changes)


@dataclass(frozen=True, slots=True, eq=False)
class ProfileAlignments:
    """Local alignments of profiles to sequence strands. ``start`` and
    ``end`` are the forward-strand half-open span of the aligned residues;
    ``score`` is in nats (:attr:`bits` in bits). Alignment ``i``'s columns
    are ``column_indptr[i]:column_indptr[i + 1]`` of ``matched`` (a residue
    fills the column, else it is deleted) and ``offset`` (along the scanned
    strand from the first aligned residue: the matched residue, or the one
    a deletion precedes); residues of the span not in a matched column are
    insertions."""

    sequence: np.ndarray
    motif: np.ndarray
    strand: np.ndarray
    start: np.ndarray
    end: np.ndarray
    score: np.ndarray
    column_indptr: np.ndarray
    matched: np.ndarray
    offset: np.ndarray

    def __len__(self) -> int:
        return len(self.score)

    @property
    def bits(self) -> np.ndarray:
        return self.score / np.log(2.0)

    def __repr__(self) -> str:
        return f"ProfileAlignments({len(self)} alignments)"


class ProfileProgram:
    """Compiled gapped profiles, aligned to every sequence strand by local
    dynamic programming. Sequences are aligned in parallel; results do not
    depend on the thread count."""

    __slots__ = ("_native", "alphabet", "background", "scoring", "strands")

    def __init__(
        self,
        native: Any,
        alphabet: Alphabet,
        background: Background,
        scoring: ProfileScoring,
        strands: Strands,
    ) -> None:
        self._native = native
        self.alphabet = alphabet
        self.background = background
        self.scoring = scoring
        self.strands = strands

    def __len__(self) -> int:
        return self._native.__len__()

    def scan(
        self,
        sequences: Sequences,
        *,
        top: int | None = None,
        names: StringVector | Seq[str] | None = None,
        limit: int | None = None,
        threads: int | None = None,
    ) -> ProfileAlignments:
        """Aligns every profile to ``sequences``.

        ``top=n`` keeps each profile's best ``n`` alignments as
        ``glam2scan`` does (a strand's next alignment is searched with the
        cells of earlier ones forbidden; ties are broken by ``names``, by
        sequence number when not given), ordered by profile then rank.
        Otherwise every alignment is kept, up to ``limit`` per strand of
        each sequence, in scan order.
        """
        if top is not None and limit is not None:
            raise ValueError("give top or limit, not both")
        s = as_strings(sequences)
        if names is not None and not isinstance(names, StringVector):
            names = StringVector.from_list(names)
        if top is not None and names is None:
            width = len(str(max(len(s) - 1, 0)))
            names = StringVector.from_list(f"{i:0{width}d}" for i in range(len(s)))
        native_names = None if names is None else (names.offsets, names.data)
        return ProfileAlignments(
            **self._native.scan(s.offsets, s.data, top, native_names, limit, threads)
        )

    def __repr__(self) -> str:
        return f"ProfileProgram({len(self)} profiles, {self.alphabet.name})"


def compile_profiles(
    motifs: MotifSource,
    background: Background,
    scoring: ProfileScoring | None = None,
    *,
    strands: Strands = "both",
    alphabet: Alphabet | None = None,
) -> ProfileProgram:
    """Compiles gapped profiles (motifs with :attr:`~omicsweave.types.Motif.gaps`)
    against the letters of ``background`` under ``scoring``
    (:meth:`ProfileScoring.glam2` by default)."""
    if strands not in ("both", "+", "-"):
        raise ValueError(f"strands must be both, + or -, not {strands!r}")
    with opened(motifs, alphabet) as (source, alphabet):
        scoring = ProfileScoring.glam2(alphabet) if scoring is None else scoring
        native = _native.ProfileProgram(
            native_source(source), alphabet, list(background.kmers), scoring, strands
        )
    return ProfileProgram(native, alphabet, background, scoring, strands)

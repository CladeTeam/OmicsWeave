"""How motif frequencies become scores and p-values."""

from __future__ import annotations

import dataclasses
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, ClassVar, Literal

import numpy as np

from .. import _native
from ..types import Alphabet
from ._background import Background

__all__ = [
    "Ambiguity",
    "AmbiguityRule",
    "ColumnPrior",
    "LogScale",
    "Quantizer",
    "Scoring",
    "Strands",
    "SubstitutionPrior",
]

Strands = Literal["both", "+", "-"]
"""The strands of each sequence a motif is tested on."""

LogScale = Literal["bits", "meme-bits", "natural", "centibits"]
"""The logarithm scores are in: ``"bits"`` (``log2``), ``"meme-bits"``
(``ln(x) · 1.44269504`` as MEME computes it), ``"natural"`` (``ln``),
``"centibits"`` (hundredths of a bit, ``ln(x) / ln(2) · 100``)."""

AmbiguityRule = Literal["skip", "mean", "sum"]
"""How a window holding an ambiguity code scores: ``"skip"`` leaves it
unscored (FIMO), ``"mean"`` scores the code by the mean of its symbols'
scores weighted by the background (MAST), ``"sum"`` by the summed
probability of its symbols (MCAST)."""


@dataclass(frozen=True, slots=True)
class Quantizer:
    """How real scores become the integers the scan and the exact null work
    with.

    - ``fimo(range)``: MEME's ``scale_score_matrix``, one offset
      ``floor(min)`` for every term and the integer scale
      ``floor(range / (max - offset))``.
    - ``precision(states)``: each term shifted by its own minimum, all
      scaled together so the best window scores ``states - 1``.
    - ``mast(range)``: MAST's ``scale_lo``, offset ``min`` and real scale
      ``range / (max - min)`` kept to 14 significant digits.
    """

    kind: Literal["fimo", "precision", "mast"]
    size: int

    @classmethod
    def fimo(cls, range: int = 1000) -> Quantizer:
        return cls("fimo", range)

    @classmethod
    def precision(cls, states: int) -> Quantizer:
        return cls("precision", states)

    @classmethod
    def mast(cls, range: int = 100) -> Quantizer:
        return cls("mast", range)


@dataclass(frozen=True, slots=True, eq=False)
class Ambiguity:
    """Ambiguity codes: each code (one character) stands for a set of
    symbols and scores by ``rule``. ``unknown`` reads every other byte
    outside the alphabet as that code (MCAST reads them as ``N``)."""

    rule: AmbiguityRule = "skip"
    codes: Mapping[str, str] = field(default_factory=dict)
    unknown: str | None = None

    def __post_init__(self) -> None:
        if self.rule not in ("skip", "mean", "sum"):
            raise ValueError(f"ambiguity rule must be skip, mean or sum, not {self.rule!r}")
        object.__setattr__(self, "codes", MappingProxyType(dict(self.codes)))

    @classmethod
    def skip(cls) -> Ambiguity:
        """No codes: a window with any symbol outside the alphabet is not scored."""
        return cls()

    @classmethod
    def iupac(cls, alphabet: Alphabet, rule: AmbiguityRule) -> Ambiguity:
        """The IUPAC nucleotide codes (both cases) of a DNA or RNA alphabet."""
        rule_, codes, unknown = _native.ambiguity_iupac(alphabet, rule)
        return cls(rule_, dict(codes), unknown)

    def with_alias(self, code: str, of: str) -> Ambiguity:
        """``code`` as another spelling of the code ``of``."""
        if of not in self.codes:
            raise ValueError(f"{of!r} is not an ambiguity code")
        return dataclasses.replace(self, codes={**self.codes, code: self.codes[of]})

    def with_unknown(self, of: str) -> Ambiguity:
        """Every byte that is neither a symbol nor a code reads as ``of``."""
        if of not in self.codes:
            raise ValueError(f"{of!r} is not an ambiguity code")
        return dataclasses.replace(self, unknown=of)

    def __repr__(self) -> str:
        unknown = f", unknown={self.unknown!r}" if self.unknown else ""
        return f"Ambiguity({self.rule!r}, {len(self.codes)} codes{unknown})"


@dataclass(frozen=True, slots=True, eq=False)
class SubstitutionPrior:
    """Smoothing through a substitution matrix: a column's frequencies
    ``f`` become ``(alpha · f + beta · g) / (alpha + beta)`` with
    ``g_i = Σ_j (f_j / b_j) · target[i, j]`` under the background letters
    ``b``. ``target`` holds the matrix's joint target frequencies
    (``alphabet × alphabet``). Columns over a reduced support are kept as
    observed."""

    kind: ClassVar[str] = "substitution"
    alpha: float
    beta: float
    target: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "target", np.ascontiguousarray(self.target, dtype=np.float64))

    @classmethod
    def mcast(cls, letters: np.ndarray) -> SubstitutionPrior:
        """MCAST's: DNA PAM-1 target frequencies at the letters ``ACGT``,
        ``alpha = 20``, ``beta = 1``."""
        alpha, beta, target = _native.substitution_mcast(
            np.ascontiguousarray(letters, dtype=np.float64)
        )
        return cls(alpha, beta, target)


ColumnPrior = SubstitutionPrior | Callable[[np.ndarray, float, np.ndarray], np.ndarray]
"""Smooths one position's frequencies in place of the pseudocount: a
:class:`SubstitutionPrior`, or a function of the column, the motif's sites
and the background letters returning the smoothed probabilities."""


@dataclass(frozen=True, slots=True, eq=False)
class Scoring:
    """Every rule that turns frequencies into scores and p-values.

    Motif rows within ``normalize_tolerance`` of 1 are kept as given;
    ``pseudocount`` (or ``column_prior``) is spread by the background over
    ``sites`` observations (``default_sites`` when the motif has none,
    truncated to an integer with ``truncate_sites``). Scores are log-odds in
    ``log``, made integer by ``quantizer``; p-values come from the exact
    null distribution of those integers (``null_budget`` bounds its
    states), under ``null_letters`` when given instead of the background.
    ``shared_reverse_null`` gives the reverse strand the forward strand's
    table. ``rounding`` keeps reported scores and p-values to that many
    significant digits (``None``: exact). ``background`` is ``"meme"`` to
    prepare it as MEME does (0-order, symmetrized on both strands, +5e-7).

    Start from a preset (:meth:`fimo`, :meth:`exact`, :meth:`mast`,
    :meth:`mcast`) and change fields with :meth:`replace`.
    """

    pseudocount: float
    default_sites: float
    truncate_sites: bool
    normalize_tolerance: float
    log: LogScale
    quantizer: Quantizer
    rounding: int | None
    background: Literal["as-given", "meme"]
    null_budget: int
    null_letters: np.ndarray | None = None
    shared_reverse_null: bool = False
    column_prior: ColumnPrior | None = None
    ambiguity: Ambiguity = field(default_factory=Ambiguity)

    def __post_init__(self) -> None:
        if self.null_letters is not None:
            object.__setattr__(
                self, "null_letters", np.ascontiguousarray(self.null_letters, dtype=np.float64)
            )

    @classmethod
    def fimo(cls) -> Scoring:
        """FIMO's: pseudocount 0.1 over 20 sites, MEME bits, 1000 levels,
        10 significant digits, MEME's background preparation."""
        return cls._from_native(_native.scoring_preset("fimo"))

    @classmethod
    def exact(cls) -> Scoring:
        """Probability semantics without MEME's numerics: exact bits, the
        background as given (any order), fine quantization, no rounding."""
        return cls._from_native(_native.scoring_preset("exact"))

    @classmethod
    def mast(
        cls,
        alphabet: Alphabet,
        *,
        letters: np.ndarray | None = None,
        mode: Literal["combine", "separate", "forward"] = "combine",
    ) -> Scoring:
        """MAST's: centibits, MAST's quantizer, p-values under ``letters``
        (MAST's built-in frequencies by default) and ambiguity codes scored
        by their mean. ``mode`` is the strand mode of the run."""
        given = None if letters is None else np.ascontiguousarray(letters, dtype=np.float64)
        return cls._from_native(_native.scoring_mast(alphabet, given, mode))

    @classmethod
    def mcast(cls, background: Background, alphabet: Alphabet) -> Scoring:
        """MCAST's: no pseudocount but a PAM substitution prior at the
        background's letters, MEME bits, ambiguity codes scored by their
        summed probability, unknown bytes read as ``N``."""
        return cls._from_native(_native.scoring_mcast(list(background.kmers), alphabet))

    @classmethod
    def _from_native(cls, d: dict[str, Any]) -> Scoring:
        kind, size = d["quantizer"]
        rule, codes, unknown = d["ambiguity"]
        prior = d["column_prior"]
        return cls(
            pseudocount=d["pseudocount"],
            default_sites=d["default_sites"],
            truncate_sites=d["truncate_sites"],
            normalize_tolerance=d["normalize_tolerance"],
            log=d["log"],
            quantizer=Quantizer(kind, size),
            rounding=d["rounding"],
            background=d["background"],
            null_budget=d["null_budget"],
            null_letters=d["null_letters"],
            shared_reverse_null=d["shared_reverse_null"],
            column_prior=None if prior is None else SubstitutionPrior(*prior),
            ambiguity=Ambiguity(rule, dict(codes), unknown),
        )

    def replace(self, **changes: Any) -> Scoring:
        """A copy with ``changes`` applied."""
        return dataclasses.replace(self, **changes)

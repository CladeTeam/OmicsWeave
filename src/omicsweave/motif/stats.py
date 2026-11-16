"""Statistics over window p-values, as replaceable parts.

A scan reports exact single-window statistics. How they become sequence or
collection statistics is a choice, made of steps each given as an object:

- :class:`BestOf`: the p-value of a sequence's best window among ``n``
  (:class:`Sidak`, :class:`MemeBestOf`)
- :class:`Combine`: one p-value from several motifs' (:class:`Product`,
  Fisher's method on the product; :class:`MinBonferroni`)
- :class:`Adjust`: E-values or q-values of a collection of p-values
  (:class:`EValue`, :class:`BenjaminiHochberg`, :class:`MemeQValue`)
- :class:`TailFit`: a right tail fitted to null samples
  (:class:`EmpiricalTail`)
- MCAST's cluster statistics, :class:`ClusterStatistics`: a null sample
  source (:class:`MemeSynthetic`, :class:`BackgroundNull`), a tail fit
  over cluster excess and GC content (:class:`GcBinnedExponential`, or
  :class:`Unconditional` around any :class:`TailFit`), and the E- and
  q-value adjustments.

The objects here run natively and can be called on numpy arrays (``q =
BenjaminiHochberg()(hits.pvalue)``). A plain Python function of the
step's shape works wherever an object is taken; it is called back from the
native layer (per value for :class:`BestOf`, per sequence for
:class:`Combine`, once per collection for :class:`Adjust`). MEME Suite's
formulas are one implementation each, not the framework's.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, ClassVar, Protocol

import numpy as np

from .. import _native
from ._defaults import default

if TYPE_CHECKING:
    from ._background import Background

__all__ = [
    "Adjust",
    "BackgroundNull",
    "BenjaminiHochberg",
    "BestOf",
    "ClusterNull",
    "ClusterStatistics",
    "ClusterTail",
    "Combine",
    "EValue",
    "EmpiricalTail",
    "GcBinnedExponential",
    "MemeBestOf",
    "MemeQValue",
    "MemeSynthetic",
    "MinBonferroni",
    "Product",
    "Sidak",
    "TailFit",
    "Unconditional",
]


class BestOf(Protocol):
    """``(p, n) -> p'``: the p-value of the best of ``n`` independent
    windows that each have p-value ``p``."""

    def __call__(self, p: float, n: float, /) -> float: ...


class Combine(Protocol):
    """``pvalues -> p``: one p-value from independent ones."""

    def __call__(self, pvalues: np.ndarray, /) -> float: ...


class Adjust(Protocol):
    """``pvalues -> values``: adjusted values (E, q, ...) in the same order."""

    def __call__(self, pvalues: np.ndarray, /) -> np.ndarray: ...


class TailFit(Protocol):
    """``samples -> tail``: a function mapping scores to the probability that
    a null score is at least as large."""

    def __call__(self, samples: np.ndarray, /) -> Callable[[float], float]: ...


class _Native:
    """A step implemented natively, read attribute by attribute."""

    __slots__ = ()
    kind: ClassVar[str]


class _NativeBestOf(_Native):
    __slots__ = ()

    def __call__(self, p: float | np.ndarray, n: float | np.ndarray) -> np.ndarray:
        p, n = np.broadcast_arrays(np.asarray(p, dtype=np.float64), np.asarray(n, dtype=np.float64))
        out = _native.stats_best_of(
            self, np.ascontiguousarray(p.ravel()), np.ascontiguousarray(n.ravel())
        )
        return out.reshape(p.shape)


class _NativeCombine(_Native):
    __slots__ = ()

    def __call__(self, pvalues: np.ndarray) -> float:
        return _native.stats_combine(self, np.ascontiguousarray(pvalues, dtype=np.float64))


class _NativeAdjust(_Native):
    __slots__ = ()

    def __call__(self, pvalues: np.ndarray) -> np.ndarray:
        return _native.stats_adjust(self, np.ascontiguousarray(pvalues, dtype=np.float64))


@dataclass(frozen=True, slots=True)
class Sidak(_NativeBestOf):
    """``1 - (1 - p) ** n``, computed without cancellation."""

    kind: ClassVar[str] = "sidak"


@dataclass(frozen=True, slots=True)
class MemeBestOf(_NativeBestOf):
    """MEME's form: ``p · n`` while that is at most 0.001, ``1 - (1 - p) ** n``
    otherwise."""

    kind: ClassVar[str] = "meme"


@dataclass(frozen=True, slots=True)
class Product(_NativeCombine):
    """Fisher's method on the product ``x`` of ``k`` p-values:
    ``x · Σ_{i<k} (-ln x)^i / i!`` (MEME's ``qfast``)."""

    kind: ClassVar[str] = "product"


@dataclass(frozen=True, slots=True)
class MinBonferroni(_NativeCombine):
    """The smallest p-value times their number, at most 1."""

    kind: ClassVar[str] = "min-bonferroni"


@dataclass(frozen=True, slots=True)
class EValue(_NativeAdjust):
    """``p · tests``; ``tests=None`` counts the p-values given."""

    kind: ClassVar[str] = "evalue"
    tests: float | None = None


@dataclass(frozen=True, slots=True)
class BenjaminiHochberg(_NativeAdjust):
    """Benjamini–Hochberg q-values, ``min_{j >= i} p_(j) · m / j`` at most 1;
    ``tests`` is ``m`` when more tests were made than p-values given."""

    kind: ClassVar[str] = "benjamini-hochberg"
    tests: float | None = None


@dataclass(frozen=True, slots=True)
class MemeQValue(_NativeAdjust):
    """MEME's q-values: Benjamini–Hochberg scaled by Storey's ``π₀``, chosen
    among ``lambdas`` values up to ``max_lambda`` by ``bootstraps``
    bootstrap rounds of ``samples`` draws (MT19937 seeded with ``seed``);
    ``π₀ = 1`` with fewer than ``min_values`` p-values."""

    kind: ClassVar[str] = "meme-qvalue"
    bootstraps: int = field(default_factory=default("meme_qvalue", "bootstraps"))
    samples: int = field(default_factory=default("meme_qvalue", "samples"))
    lambdas: int = field(default_factory=default("meme_qvalue", "lambdas"))
    max_lambda: float = field(default_factory=default("meme_qvalue", "max_lambda"))
    min_values: int = field(default_factory=default("meme_qvalue", "min_values"))
    seed: int = field(default_factory=default("meme_qvalue", "seed"))


@dataclass(frozen=True, slots=True)
class EmpiricalTail(_Native):
    """The empirical tail, ``(1 + #{samples >= x}) / (1 + n)``."""

    kind: ClassVar[str] = "empirical"

    def __call__(self, samples: np.ndarray) -> Callable[[float], float]:
        sorted_ = np.sort(np.asarray(samples, dtype=np.float64))
        if not len(sorted_):
            raise ValueError("no samples to fit a tail to")
        n = len(sorted_)

        def tail(x: float) -> float:
            return (1 + n - int(np.searchsorted(sorted_, x, side="left"))) / (1 + n)

        return tail


# --- MCAST's cluster statistics -------------------------------------------


class ClusterNull(_Native):
    """A source of null cluster scores."""

    __slots__ = ()


class ClusterTail(_Native):
    """A null fitted to cluster samples: p-value of a score excess at a GC
    content."""

    __slots__ = ()


@dataclass(frozen=True, slots=True)
class MemeSynthetic(ClusterNull):
    """MCAST's null: rounds of one random ``round_length`` sequence each, at
    a GC content drawn uniformly between the data's clusters' lowest and
    highest, until ``min_rounds`` rounds found ``want`` clusters; it gives
    up after ``give_up_after`` rounds when far short and stops at
    ``max_rounds``. ``max_samples`` clusters are kept by reservoir
    sampling. One MT19937 stream seeded with ``seed`` draws in MCAST's
    order. Needs an ``ACGT`` alphabet."""

    kind: ClassVar[str] = "meme-synthetic"
    seed: int = 0
    round_length: int = field(default_factory=default("meme_synthetic", "round_length"))
    min_rounds: int = field(default_factory=default("meme_synthetic", "min_rounds"))
    max_rounds: int = field(default_factory=default("meme_synthetic", "max_rounds"))
    give_up_after: int = field(default_factory=default("meme_synthetic", "give_up_after"))
    want: float = field(default_factory=default("meme_synthetic", "want"))
    max_samples: int = field(default_factory=default("meme_synthetic", "max_samples"))


@dataclass(frozen=True, slots=True, eq=False)
class BackgroundNull(ClusterNull):
    """Null clusters from ``sequences`` sequences of ``length`` letters drawn
    from ``background``."""

    kind: ClassVar[str] = "background"
    background: Background
    length: int
    sequences: int
    seed: int = 0


@dataclass(frozen=True, slots=True)
class GcBinnedExponential(ClusterTail):
    """MCAST's tail: exponential tails ``exp(-x / μ)`` fitted in GC-content
    bins of at least ``min_scores`` samples (a bin closes once its running
    mean GC moves by the ratio ``max_ratio``), ``μ`` interpolated between
    bins. Samples with fewer than ``2 · max_gap`` letters around them in
    their segment are dropped."""

    kind: ClassVar[str] = "gc-binned-exponential"
    max_gap: int
    min_scores: int = field(default_factory=default("gc_binned_exponential", "min_scores"))
    max_ratio: float = field(default_factory=default("gc_binned_exponential", "max_ratio"))


@dataclass(frozen=True, slots=True)
class Unconditional(ClusterTail):
    """``fit`` applied to the excess scores alone, ignoring GC content."""

    kind: ClassVar[str] = "unconditional"
    fit: TailFit | EmpiricalTail = field(default_factory=EmpiricalTail)


@dataclass(frozen=True, slots=True, eq=False)
class ClusterStatistics:
    """How cluster p-, E- and q-values are computed: null samples, a tail
    fitted to them, and the adjustments of the resulting p-values."""

    null: ClusterNull
    tail: ClusterTail
    evalue: Adjust | EValue = field(default_factory=EValue)
    qvalue: Adjust | MemeQValue = field(default_factory=MemeQValue)

    @classmethod
    def mcast(cls, *, seed: int = 0, max_gap: int = 50) -> ClusterStatistics:
        """MCAST's: synthetic null with ``--seed``, GC-binned exponential
        tails, E-values over the clusters and MEME's q-values."""
        return cls(MemeSynthetic(seed), GcBinnedExponential(max_gap), EValue(), MemeQValue())

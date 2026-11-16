"""What a scan keeps (reductions) and what it returns (results).

A reduction names one product of a scan and the level it reduces to:
windows (:class:`Hits`, :class:`Scores`), one sequence (:class:`Best`,
:class:`Counts`, :class:`Clusters`), or the collection (:class:`Mast`).
Any number of them share one pass over the sequences. Results are columns
of numpy arrays; ``strand`` columns are ``int8``: ``1`` forward, ``-1``
reverse (``0`` where strands are combined).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, ClassVar, Literal

import numpy as np

from ._defaults import default
from .stats import (
    Adjust,
    BestOf,
    ClusterStatistics,
    Combine,
    EValue,
    MemeBestOf,
    Product,
)

__all__ = [
    "Best",
    "BestTable",
    "ClusterTable",
    "Clusters",
    "CountMatrix",
    "Counts",
    "HitTable",
    "Hits",
    "Mast",
    "MastHits",
    "MastMotifs",
    "MastReport",
    "Reduction",
    "ScoreTracks",
    "Scores",
]

StrandMode = Literal["combine", "separate", "forward"]


def _threshold(pvalue: float | None, score: float | None) -> tuple[str, float]:
    if pvalue is not None and score is not None:
        raise ValueError("give a p-value or a score threshold, not both")
    if score is not None:
        return ("score", float(score))
    return ("pvalue", 1e-4 if pvalue is None else float(pvalue))


class Reduction:
    """A product of a scan. ``part``, where a reduction has it, limits it to
    the motifs of one part of a combined program
    (:func:`~omicsweave.motif.combine`)."""

    __slots__ = ()
    kind: ClassVar[str]


@dataclass(frozen=True, slots=True)
class Hits(Reduction):
    """Every window passing a threshold: a p-value at most ``pvalue`` (1e-4
    when neither is given) or a score at least ``score``. Returns a
    :class:`HitTable`."""

    kind: ClassVar[str] = "hits"
    pvalue: float | None = None
    score: float | None = None
    part: int | None = None

    @property
    def threshold(self) -> tuple[str, float]:
        return _threshold(self.pvalue, self.score)


@dataclass(frozen=True, slots=True)
class Best(Reduction):
    """Each sequence's best window per motif (over the strands scanned).
    Returns a :class:`BestTable`."""

    kind: ClassVar[str] = "best"
    part: int | None = None


@dataclass(frozen=True, slots=True)
class Counts(Reduction):
    """Windows passing a threshold (as :class:`Hits`), counted per sequence
    and motif over both strands. Returns a :class:`CountMatrix`."""

    kind: ClassVar[str] = "counts"
    pvalue: float | None = None
    score: float | None = None
    part: int | None = None

    @property
    def threshold(self) -> tuple[str, float]:
        return _threshold(self.pvalue, self.score)


@dataclass(frozen=True, slots=True)
class Scores(Reduction):
    """Every window's score of every motif, densely. Returns
    :class:`ScoreTracks`."""

    kind: ClassVar[str] = "scores"


@dataclass(frozen=True, slots=True, eq=False)
class Mast(Reduction):
    """MAST's sequence scores: per sequence (per strand with
    ``mode="separate"``), each motif's best window becomes a sequence
    p-value through ``best_of``, the motifs' p-values combine through
    ``combine`` and E-values come from ``adjust`` over every sequence.
    Windows with p-value below ``hit_threshold`` tile each sequence; the
    sequences with E-value at most ``evalue_threshold`` are ranked. The
    program should be compiled with :meth:`Scoring.mast` under the same
    ``mode``. Returns a :class:`MastReport`."""

    kind: ClassVar[str] = "mast"
    mode: StrandMode = "combine"
    hit_threshold: float = field(default_factory=default("mast", "hit_threshold"))
    evalue_threshold: float = field(default_factory=default("mast", "evalue_threshold"))
    best_of: BestOf | MemeBestOf = field(default_factory=MemeBestOf)
    combine: Combine | Product = field(default_factory=Product)
    adjust: Adjust | EValue = field(default_factory=EValue)
    part: int | None = None


@dataclass(frozen=True, slots=True, eq=False)
class Clusters(Reduction):
    """MCAST's clusters: windows with p-value at most ``pthresh`` are hits,
    scored ``-log2(p / pthresh)``; a gap costs ``threshold / max_gap`` per
    letter, and a cluster is a maximal run of hits whose score less its
    gaps exceeds the threshold ``D`` (from ``pthresh``, ``max_gap`` and
    ``egcost`` unless given). Long sequences are solved in segments of
    ``segment`` letters overlapping by ``overlap``. With ``statistics``,
    clusters get p-, E- and q-values. The program should be compiled with
    :meth:`Scoring.mcast` on both strands. Returns a :class:`ClusterTable`."""

    kind: ClassVar[str] = "clusters"
    pthresh: float = field(default_factory=default("clusters", "pthresh"))
    max_gap: int = field(default_factory=default("clusters", "max_gap"))
    egcost: float = field(default_factory=default("clusters", "egcost"))
    threshold: float | None = None
    segment: int | None = None
    overlap: int = field(default_factory=default("clusters", "overlap"))
    statistics: ClusterStatistics | None = None
    part: int | None = None


# --- results ---------------------------------------------------------------


@dataclass(frozen=True, slots=True, eq=False)
class HitTable:
    """Windows, sorted by motif, sequence, start and strand. ``start`` and
    ``end`` are 0-based half-open on the forward strand; ``score`` is the
    reported (rounded) score and ``pvalue`` its exact p-value."""

    sequence: np.ndarray
    motif: np.ndarray
    start: np.ndarray
    end: np.ndarray
    strand: np.ndarray
    score: np.ndarray
    pvalue: np.ndarray

    def __len__(self) -> int:
        return len(self.score)

    def __repr__(self) -> str:
        return f"HitTable({len(self)} hits)"


@dataclass(frozen=True, slots=True, eq=False)
class BestTable:
    """Per (sequence, motif) that scored, the best window: higher score
    first, then smaller start, then forward strand. Sorted by sequence and
    motif."""

    sequence: np.ndarray
    motif: np.ndarray
    start: np.ndarray
    strand: np.ndarray
    score: np.ndarray
    pvalue: np.ndarray

    def __len__(self) -> int:
        return len(self.score)

    def __repr__(self) -> str:
        return f"BestTable({len(self)} rows)"


@dataclass(frozen=True, slots=True, eq=False)
class CountMatrix:
    """A sparse sequence × motif matrix in CSR layout: sequence ``i``'s
    nonzero counts are ``counts[indptr[i]:indptr[i + 1]]`` at motifs
    ``motifs[...]``."""

    indptr: np.ndarray
    motifs: np.ndarray
    counts: np.ndarray
    shape: tuple[int, int]

    def __repr__(self) -> str:
        return f"CountMatrix({self.shape[0]} × {self.shape[1]}, {len(self.counts)} nonzero)"

    def to_dense(self) -> np.ndarray:
        out = np.zeros(self.shape, dtype=np.uint32)
        rows = np.repeat(np.arange(self.shape[0]), np.diff(self.indptr).astype(np.int64))
        out[rows, self.motifs] = self.counts
        return out

    def to_scipy(self) -> Any:
        """A ``scipy.sparse.csr_array`` over the same buffers."""
        from scipy.sparse import csr_array

        return csr_array((self.counts, self.motifs, self.indptr), shape=self.shape)


@dataclass(frozen=True, slots=True, eq=False)
class ScoreTracks:
    """Every window's reported score (``float32``, NaN where no window
    scores). Sequence ``i`` is ``values[offsets[i]:offsets[i + 1]]``, a
    ``(motifs, strands, length_i)`` array indexed by its windows' starts;
    strands in program order, forward first."""

    offsets: np.ndarray
    values: np.ndarray
    motifs: int
    strands: int
    lengths: np.ndarray | None = None
    """Sequence lengths; retained even when the program has no motifs."""

    def __len__(self) -> int:
        return len(self.offsets) - 1

    def __repr__(self) -> str:
        return f"ScoreTracks({len(self)} sequences, {self.motifs} motifs, {self.strands} strands)"

    def __getitem__(self, i: int) -> np.ndarray:
        n = len(self)
        if i < 0:
            i += n
        if not 0 <= i < n:
            raise IndexError("sequence index out of range")
        block = self.values[self.offsets[i] : self.offsets[i + 1]]
        size = self.motifs * self.strands
        length = (
            int(self.lengths[i]) if self.lengths is not None else len(block) // size if size else 0
        )
        return block.reshape(self.motifs, self.strands, length)


@dataclass(frozen=True, slots=True, eq=False)
class MastMotifs:
    """Per report row ``r``, rows ``indptr[r]:indptr[r + 1]``: each motif's
    best window score, its position p-value and the sequence p-value
    ``best_of`` made of it."""

    indptr: np.ndarray
    motif: np.ndarray
    score: np.ndarray
    pvalue: np.ndarray
    sequence_pvalue: np.ndarray


@dataclass(frozen=True, slots=True, eq=False)
class MastHits:
    """Per report row ``r``, rows ``indptr[r]:indptr[r + 1]``: the
    non-overlapping strong windows tiling it, by start."""

    indptr: np.ndarray
    motif: np.ndarray
    start: np.ndarray
    end: np.ndarray
    strand: np.ndarray
    score: np.ndarray
    pvalue: np.ndarray


@dataclass(frozen=True, slots=True, eq=False)
class MastReport:
    """One row per sequence (per strand of it with ``mode="separate"``, in
    input order): its combined ``pvalue`` and ``evalue``, its motifs'
    ``best`` and its tiling ``hits``. ``ranked`` lists the rows within the
    E-value threshold, by p-value, then input order, then strand."""

    sequence: np.ndarray
    strand: np.ndarray
    pvalue: np.ndarray
    evalue: np.ndarray
    ranked: np.ndarray
    best: MastMotifs
    hits: MastHits

    def __len__(self) -> int:
        return len(self.sequence)

    def __repr__(self) -> str:
        return f"MastReport({len(self)} rows, {len(self.ranked)} ranked)"


@dataclass(frozen=True, slots=True, eq=False)
class ClusterTable:
    """Clusters: by increasing p-value with statistics (then decreasing
    score, decreasing GC, discovery order), in discovery order without.
    ``id`` is the 1-based discovery rank; ``start`` / ``end`` 0-based
    half-open; ``excess`` the score less the threshold ``D``; ``hits`` the
    windows inside; ``gc`` the GC fraction of the cluster and 500 letters
    either side. ``pvalue``, ``evalue`` and ``qvalue`` are NaN when
    ``statistics`` is false."""

    id: np.ndarray
    sequence: np.ndarray
    start: np.ndarray
    end: np.ndarray
    score: np.ndarray
    excess: np.ndarray
    hits: np.ndarray
    gc: np.ndarray
    segment_start: np.ndarray
    segment_len: np.ndarray
    pvalue: np.ndarray
    evalue: np.ndarray
    qvalue: np.ndarray
    statistics: bool
    threshold: float
    """The cluster threshold ``D`` used."""

    def __len__(self) -> int:
        return len(self.id)

    def __repr__(self) -> str:
        return f"ClusterTable({len(self)} clusters, statistics={self.statistics})"


def build_result(d: dict[str, Any]) -> Any:
    """A native result (a dict naming its kind) as its result class."""
    kind = d.pop("kind")
    match kind:
        case "hits":
            return HitTable(**d)
        case "best":
            return BestTable(**d)
        case "counts":
            return CountMatrix(**d)
        case "scores":
            return ScoreTracks(**d)
        case "mast":
            best = MastMotifs(**d.pop("best"))
            hits = MastHits(**d.pop("hits"))
            return MastReport(best=best, hits=hits, **d)
        case "clusters":
            return ClusterTable(**d)
    raise ValueError(f"unknown result kind {kind!r}")

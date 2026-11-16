"""What the columns of a count matrix are."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

from ..formats import PathLike

if TYPE_CHECKING:
    from ..annotation import IntervalSet

__all__ = [
    "Bins",
    "Features",
    "GeneActivity",
    "GeneExons",
    "GeneSpans",
    "Peaks",
    "Promoters",
    "Regions",
]


@dataclass(frozen=True, slots=True)
class GeneExons:
    """Genes by their merged exons (featureCounts meta-features, STARsolo
    ``Gene``); needs the task's annotation."""

    kind: ClassVar[str] = "gene-exons"


@dataclass(frozen=True, slots=True)
class GeneSpans:
    """Genes by their whole span, introns included (STARsolo ``GeneFull``);
    needs the task's annotation."""

    kind: ClassVar[str] = "gene-spans"


@dataclass(frozen=True, slots=True)
class Promoters:
    """A window around each gene's TSS, ``upstream`` and ``downstream``
    bases on the gene's strand, clipped to the chromosomes."""

    upstream: int = 2000
    downstream: int = 200
    kind: ClassVar[str] = "promoters"


@dataclass(frozen=True, slots=True)
class GeneActivity:
    """Each gene body extended ``upstream`` and ``downstream`` (Signac gene
    activity), clipped to the chromosomes."""

    upstream: int = 2000
    downstream: int = 0
    kind: ClassVar[str] = "gene-activity"


@dataclass(frozen=True, slots=True)
class Peaks:
    """Every feature of a bigBed file, one column each, named
    ``chrom:start-end`` (BED converts with ``ow convert peaks.bed
    peaks.bb``)."""

    path: PathLike
    kind: ClassVar[str] = "peaks"


@dataclass(frozen=True, slots=True)
class Bins:
    """Tiles of ``size`` bases over every reference of the first input."""

    size: int
    kind: ClassVar[str] = "bins"


@dataclass(frozen=True, slots=True, eq=False)
class Regions:
    """Every row of an :class:`~omicsweave.annotation.IntervalSet` as its
    own column (introns, exon bins, promoters of a subset, or
    :meth:`IntervalSet.from_arrays <omicsweave.annotation.IntervalSet.from_arrays>`),
    with id ``chrom:start-end`` and the row's name; ``feature_type`` is the
    third column of ``features.tsv.gz``."""

    intervals: IntervalSet
    feature_type: str = "Regions"
    kind: ClassVar[str] = "regions"


Features = GeneExons | GeneSpans | Promoters | GeneActivity | Peaks | Bins | Regions

_BY_KIND: dict[str, type] = {
    c.kind: c for c in (GeneExons, GeneSpans, Promoters, GeneActivity, Peaks, Bins)
}


def from_native(values: dict[str, object]) -> Features:
    """A preset's features: a dict naming its ``kind``."""
    values = dict(values)
    return _BY_KIND[values.pop("kind")](**values)  # type: ignore[return-value]

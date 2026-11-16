"""Signal by locus, from bigWig: per position, in bins, or as stored runs."""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from types import MappingProxyType

import numpy as np

from .. import _native
from ..cache import BlockCache, BlockKey
from ..cache import native as native_cache
from ..errors import InvalidError, UnknownChromError
from ..formats import BigWig, BigWigSummary, PathLike
from ..types import LocusLike, RegionLike, Signal, SignalBins, SignalIntervals, Stat
from ..types._genome import Locus, as_locus, as_region, native_key
from ._base import CachedReader, Reader, resolve

__all__ = ["StrandedTrackReader", "TrackReader"]


def _stats(stats: Stat | str | Iterable[Stat | str]) -> tuple[Stat, ...]:
    if isinstance(stats, str):
        return (Stat(stats),)
    out = tuple(Stat(s) for s in stats)
    if not out:
        raise InvalidError("bins: no statistic asked for", what="bins", detail="no statistic")
    return out


def _bin_count(length: int, resolution: int, phase: int) -> int:
    if not isinstance(resolution, (int, np.integer)) or resolution <= 0:
        raise InvalidError(
            "bins: resolution must be a positive integer", what="bins", detail="resolution"
        )
    if not 0 <= phase < resolution:
        detail = f"phase {phase} must lie in [0, {resolution})"
        raise InvalidError(f"bins: {detail}", what="bins", detail=detail)
    return 0 if length == 0 else -(-(length + phase) // resolution)


class TrackReader(CachedReader):
    """Reads the signal of a bigWig file by locus: a
    :class:`~omicsweave.types.Region` or :class:`~omicsweave.types.Spliced`
    spans, 5' to 3' of the key's strand."""

    def __init__(self, source: BigWig | PathLike, *, cache: BlockCache | None = None) -> None:
        file = resolve(source, (BigWig,), "TrackReader")
        native = _native.TrackReader(os.fspath(file), native_cache(cache))
        self._open(native, file)
        self.source = file
        self._chroms = MappingProxyType(dict(native.chroms()))
        summary = native.summary()
        self._summary = None if summary is None else BigWigSummary(*summary)

    def set_cache(self, cache: BlockCache | None) -> None:
        """Attach external decoded blocks, or detach them with None."""
        self._native.set_cache(native_cache(cache))

    def describe_cache(self, locus: LocusLike) -> tuple[BlockKey, ...]:
        """The full-resolution blocks this locus depends on, each once."""
        return tuple(self._native.describe_cache(native_key(as_locus(locus))))

    def prefetch_cache(self, locus: LocusLike) -> None:
        """Prepare blocks without constructing the query's final output."""
        self._native.prefetch_cache(native_key(as_locus(locus)))

    @property
    def chroms(self) -> Mapping[str, int]:
        """Chromosome lengths by name."""
        return self._chroms

    @property
    def summary(self) -> BigWigSummary | None:
        return self._summary

    @property
    def zoom_reductions(self) -> tuple[int, ...]:
        """Bases per record of each stored zoom level, finest first."""
        return tuple(self._native.zoom_reductions())

    def signal(self, locus: LocusLike, out: np.ndarray | None = None) -> Signal:
        """One value per position of ``locus``; it may run past its
        chromosome's end, where there is no data (NaN).

        By default, returns a new array independent of later queries and
        the reader's lifetime. Optional ``out`` (writable, contiguous
        ``float32``, one element per position) receives the values in
        place; the returned signal references that same array.
        """
        key = as_locus(locus)
        if out is None:
            if self.closed:
                raise ValueError("I/O operation on closed file")
            if key.chrom not in self._chroms:
                raise UnknownChromError(f"unknown chromosome {key.chrom!r}", chrom=key.chrom)
            out = np.empty(len(key), dtype=np.float32)
        self._native.values(native_key(key), out)
        return Signal(key, out)

    def bins(
        self,
        locus: LocusLike,
        resolution: int,
        stats: Stat | str | Iterable[Stat | str] = Stat.MEAN,
        *,
        phase: int = 0,
        exact: bool = True,
        out: np.ndarray | None = None,
    ) -> SignalBins:
        """The ``stats`` of ``locus`` in bins of ``resolution`` positions,
        bins 5' to 3' (see :class:`~omicsweave.types.SignalBins`).

        All statistics come from one pass over the data. ``exact`` summarizes
        the full-resolution data, bit for bit as libBigWig's
        ``bwStatsFromFull`` on full bins; ``exact=False`` reads the coarsest
        zoom level finer than half a bin, as ``bwStats`` picks it, weighting
        each zoom record by the fraction of it inside the bin: an estimate
        that reads far less for wide bins (and the full data when no level
        qualifies). Results are allocated automatically and remain valid
        independently of later queries and the reader's lifetime. Optional
        ``out`` (writable, C-contiguous ``float32``, ``(len(stats), bins)``)
        receives the values in place; the result references that same array.
        """
        key = as_locus(locus)
        stats = _stats(stats)
        n = _bin_count(len(key), resolution, phase)
        if out is None:
            if self.closed:
                raise ValueError("I/O operation on closed file")
            if key.chrom not in self._chroms:
                raise UnknownChromError(f"unknown chromosome {key.chrom!r}", chrom=key.chrom)
            out = np.empty((len(stats), n), dtype=np.float32)
        self._native.bins(native_key(key), resolution, phase, [str(s) for s in stats], exact, out)
        return SignalBins(key, int(resolution), int(phase), stats, out)

    def intervals(self, region: RegionLike) -> SignalIntervals:
        """The runs of constant signal in ``region``, as stored, in genomic
        order whatever the strand."""
        r = as_region(region)
        return SignalIntervals(r, *self._native.intervals(r.chrom, r.start, r.end))


class StrandedTrackReader(Reader):
    """Reads a stranded signal stored as two bigWig files, one per strand
    (RNA-seq, CAGE, PRO-seq coverage). A key on ``"+"`` reads ``plus``; a key
    on ``"-"`` reads ``minus``, 5' to 3' of the minus strand. Values come as
    stored (a minus file holding negative values keeps them). Unstranded
    keys are refused: they do not say which file to read."""

    def __init__(
        self,
        plus: BigWig | PathLike,
        minus: BigWig | PathLike,
        *,
        cache: BlockCache | None = None,
    ) -> None:
        self.plus = TrackReader(plus, cache=cache)
        self.minus = TrackReader(minus, cache=cache)
        self._open(_Pair(self.plus, self.minus), self.plus.source, self.minus.source)

    def _side(self, locus: LocusLike) -> tuple[TrackReader, Locus]:
        key = as_locus(locus)
        if key.strand == "+":
            return self.plus, key
        if key.strand == "-":
            return self.minus, key
        raise InvalidError(
            f"a stranded track reads keys on + or -, not {key.strand!r}",
            what="stranded track",
            detail="key strand",
        )

    @property
    def chroms(self) -> Mapping[str, int]:
        return self.plus.chroms

    def signal(self, locus: LocusLike, out: np.ndarray | None = None) -> Signal:
        """The values of the key's strand; see :meth:`TrackReader.signal`."""
        track, key = self._side(locus)
        return track.signal(key, out)

    def bins(
        self,
        locus: LocusLike,
        resolution: int,
        stats: Stat | str | Iterable[Stat | str] = Stat.MEAN,
        *,
        phase: int = 0,
        exact: bool = True,
        out: np.ndarray | None = None,
    ) -> SignalBins:
        """The bins of the key's strand; see :meth:`TrackReader.bins`."""
        track, key = self._side(locus)
        return track.bins(key, resolution, stats, phase=phase, exact=exact, out=out)

    def intervals(self, region: RegionLike) -> SignalIntervals:
        """The stored runs of the region's strand, in genomic order."""
        track, key = self._side(as_region(region))
        return track.intervals(key)  # type: ignore[arg-type]


class _Pair:
    """Closes both strands' readers."""

    def __init__(self, plus: TrackReader, minus: TrackReader) -> None:
        self._readers = (plus, minus)

    def close(self) -> None:
        for r in self._readers:
            r.close()

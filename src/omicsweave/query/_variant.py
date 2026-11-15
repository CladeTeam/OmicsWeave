"""Genotypes by region from BCF."""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from types import MappingProxyType

from .. import _native
from ..cache import BlockCache, BlockKey
from ..cache import native as native_cache
from ..formats import Bcf, Csi, PathLike, detect
from ..types import Genotypes, RegionLike, StringVector
from ..types._genome import as_region
from ._base import CachedReader, resolve

__all__ = ["VariantReader"]


class VariantReader(CachedReader):
    """Reads the records of a BCF file by region, through its CSI index, as
    htslib's ``bcf_itr_queryi`` does, and the genotypes of samples at them.

    ``index`` defaults to ``<bcf>.csi`` beside the file. A missing index
    raises :class:`~omicsweave.errors.IoError`: it is written with the BCF
    (``ow convert cohort.vcf.gz cohort.bcf --index``), never on open. VCF
    text is not read by region: convert it to BCF first.

    ``cache`` optionally supplies caller-owned decoded BGZF blocks. The BCF
    and its index must remain unchanged throughout the reader session.
    """

    def __init__(
        self,
        source: Bcf | PathLike,
        index: Csi | PathLike | None = None,
        *,
        cache: BlockCache | None = None,
    ) -> None:
        file = resolve(source, (Bcf,), "VariantReader")
        native = _native.VariantReader(
            os.fspath(file), None if index is None else os.fspath(index), native_cache(cache)
        )
        self._open(native, file, index)
        self.source = file
        self._chroms = MappingProxyType(dict(native.chroms()))
        self._samples = tuple(native.samples())
        self._sample_ids = {}
        for i, name in enumerate(self._samples):
            self._sample_ids.setdefault(name, i)
        self.index = detect(native.index_path)

    def set_cache(self, cache: BlockCache | None) -> None:
        """Attach caller-owned BGZF blocks, or detach them with None."""
        self._native.set_cache(native_cache(cache))

    def describe_cache(self, region: RegionLike) -> tuple[BlockKey, ...]:
        """The physical BGZF dependencies of the indexed chunks."""
        r = as_region(region)
        return tuple(self._native.describe_cache(r.chrom, r.start, r.end))

    def prefetch_cache(self, region: RegionLike) -> None:
        """Load blocks without constructing records or a genotype matrix."""
        r = as_region(region)
        self._native.prefetch_cache(r.chrom, r.start, r.end)

    def stream_stats(self) -> dict[str, int]:
        """Physical reads, inflations and buffer movement in this context."""
        return self._native.stream_stats()

    def reset_stream_stats(self) -> None:
        """Reset physical-work counters, retaining blocks and query state."""
        self._native.reset_stream_stats()

    @property
    def chroms(self) -> Mapping[str, int | None]:
        """Contig lengths by name, as the header declares them (``None``
        where it gives none)."""
        return self._chroms

    @property
    def samples(self) -> tuple[str, ...]:
        """Sample names in header order."""
        return self._samples

    def genotypes(
        self, region: RegionLike, samples: str | Iterable[str] | None = None
    ) -> Genotypes:
        """The records overlapping ``region`` and the genotypes of
        ``samples`` (names, in the order given; default: every sample)."""
        r = as_region(region)
        if samples is None:
            names = self._samples
            ids = None
        else:
            names = (samples,) if isinstance(samples, str) else tuple(samples)
            try:
                ids = [self._sample_ids[name] for name in names]
            except KeyError as e:
                raise KeyError(f"no sample {e.args[0]!r} in {self.source.path}") from None
        start, end, ids_, ref, alt, ploidy, gt, phased = self._native.genotypes(
            r.chrom, r.start, r.end, ids
        )
        n = len(start)
        return Genotypes(
            r,
            names,
            start,
            end,
            StringVector(*ids_),
            StringVector(*ref),
            StringVector(*alt),
            gt.reshape(n, len(names), ploidy),
            phased.reshape(n, len(names)),
        )

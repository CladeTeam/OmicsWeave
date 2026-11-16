"""Records and genotypes by region from BCF."""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from types import MappingProxyType
from typing import Any

from .. import _native
from ..cache import BlockCache
from ..cache import native as native_cache
from ..formats import Bcf, Csi, PathLike, detect
from ..types import Genotypes, RegionLike, StringVector, Variants
from ..types._genome import as_region
from ._base import CachedReader, resolve

__all__ = ["VariantReader"]


class VariantReader(CachedReader):
    """Reads the records of a BCF file by region, through its CSI index, as
    htslib's ``bcf_itr_queryi`` does: :meth:`fetch` with QUAL, FILTER and
    any INFO / FORMAT fields, optionally filtered by a bcftools expression;
    :meth:`genotypes` the genotype matrix alone, decoding nothing else.

    ``index`` defaults to the one htslib finds beside the file
    (``<bcf>.csi``, ``<name>.csi``). A missing index raises
    :class:`~omicsweave.errors.IoError`: indexes are built explicitly
    (:meth:`Bcf.build_index <omicsweave.formats.Bcf.build_index>`, or with
    the BCF: ``ow convert cohort.vcf.gz cohort.bcf --index``), never on
    open. VCF text is not read by region: convert it to BCF first.

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

    @property
    def chroms(self) -> Mapping[str, int | None]:
        """Contig lengths by name, as the header declares them (``None``
        where it gives none)."""
        return self._chroms

    @property
    def samples(self) -> tuple[str, ...]:
        """Sample names in header order."""
        return self._samples

    def _samples_of(
        self, samples: str | Iterable[str] | None
    ) -> tuple[tuple[str, ...], list[int] | None]:
        if samples is None:
            return self._samples, None
        names = (samples,) if isinstance(samples, str) else tuple(samples)
        try:
            return names, [self._sample_ids[name] for name in names]
        except KeyError as e:
            raise KeyError(f"no sample {e.args[0]!r} in {self.source.path}") from None

    def genotypes(
        self, region: RegionLike, samples: str | Iterable[str] | None = None
    ) -> Genotypes:
        """The records overlapping ``region`` and the genotypes of
        ``samples`` (names, in the order given; default: every sample),
        decoding nothing but the fixed fields and ``GT``."""
        r = as_region(region)
        names, ids = self._samples_of(samples)
        start, end, id_, ref, alt, ploidy, gt, phased = self._native.genotypes(
            r.chrom, r.start, r.end, ids
        )
        columns = dict(
            start=start, end=end, id=id_, ref=ref, alt=alt, ploidy=ploidy, gt=gt, phased=phased
        )
        return _genotypes(r, names, columns)

    def fetch(
        self,
        region: RegionLike,
        samples: str | Iterable[str] | None = None,
        *,
        info: Iterable[str] = (),
        format: Iterable[str] = (),
        include: str | None = None,
        exclude: str | None = None,
    ) -> Variants:
        """The records overlapping ``region`` with QUAL, FILTER, the
        ``info`` and ``format`` fields named (as the header declares them)
        and the genotypes of ``samples``, in one pass. ``include`` keeps the
        records a bcftools expression passes (``-i``), ``exclude`` leaves
        them out (``-e``); e.g. ``include="QUAL>30 && INFO/DP>10"``."""
        r = as_region(region)
        names, ids = self._samples_of(samples)
        info, format = _names(info, "info"), _names(format, "format")
        c = self._native.records(
            r.chrom, r.start, r.end, ids, info, format, include=include, exclude=exclude
        )
        g = _genotypes(r, names, c)
        n, k = len(g), len(names)
        return Variants(
            g,
            c["qual"],
            StringVector(*c["filters"]),
            MappingProxyType({f[0]: _field(f, (n,)) for f in c["info"]}),
            MappingProxyType({f[0]: _field(f, (n, k)) for f in c["format"]}),
        )


def _names(fields: Iterable[str], what: str) -> list[str]:
    if isinstance(fields, str):
        fields = (fields,)
    out = list(fields)
    if len(set(out)) != len(out):
        raise ValueError(f"{what}: a field is named twice")
    return out


def _genotypes(region: Any, names: tuple[str, ...], c: Mapping[str, Any]) -> Genotypes:
    n = len(c["start"])
    return Genotypes(
        region,
        names,
        c["start"],
        c["end"],
        StringVector(*c["id"]),
        StringVector(*c["ref"]),
        StringVector(*c["alt"]),
        c["gt"].reshape(n, len(names), c["ploidy"]),
        c["phased"].reshape(n, len(names)),
    )


def _field(field: tuple[str, str, int, Any], lead: tuple[int, ...]) -> Any:
    """A field's values: flags as they come, numbers shaped ``lead`` (records,
    and samples for FORMAT) by the width, strings as a vector."""
    _, kind, width, values = field
    if kind == "flag":
        return values
    if kind == "str":
        return StringVector(*values)
    return values.reshape(*lead, width)

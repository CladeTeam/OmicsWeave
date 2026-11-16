"""Opening the reader a file's format calls for."""

from __future__ import annotations

from ..errors import InvalidError
from ..formats import (
    Bam,
    Bcf,
    BgzipFasta,
    BigBed,
    BigWig,
    Cram,
    Fasta,
    File,
    Hic,
    PathLike,
    TwoBit,
    detect,
)
from ._alignment import AlignmentReader
from ._base import conversion
from ._contact import ContactReader
from ._interval import IntervalReader
from ._sequence import SequenceReader
from ._track import TrackReader
from ._variant import VariantReader

__all__ = ["open"]


def open(
    source: File | PathLike, **options: object
) -> (
    SequenceReader | TrackReader | IntervalReader | ContactReader | AlignmentReader | VariantReader
):
    """The reader of ``source``'s format (a path is detected first);
    ``options`` go to that reader.

    Formats that are not read by region (VCF, SAM, BED and signal text,
    plain gzip FASTA, cooler, pairs) raise
    :class:`~omicsweave.errors.InvalidError` naming the conversion that makes
    them readable.
    """
    file = source if isinstance(source, File) else detect(source)
    detail = conversion(file)
    if detail is not None:
        raise InvalidError(f"{file.path}: {detail}", what=str(file.path), detail=detail)
    match file:
        case Fasta() | BgzipFasta() | TwoBit():
            return SequenceReader(file, **options)  # type: ignore[arg-type]
        case BigWig():
            return TrackReader(file, **options)  # type: ignore[arg-type]
        case BigBed():
            return IntervalReader(file, **options)  # type: ignore[arg-type]
        case Hic():
            return ContactReader(file, **options)  # type: ignore[arg-type]
        case Bam() | Cram():
            return AlignmentReader(file, **options)  # type: ignore[arg-type]
        case Bcf():
            return VariantReader(file, **options)  # type: ignore[arg-type]
    raise TypeError(f"no reader reads {type(file).__name__}")

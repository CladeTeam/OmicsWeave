"""Opening the reader a file's format calls for."""

from __future__ import annotations

from ..errors import InvalidError
from ..formats import (
    Bcf,
    Bed,
    BgzipFasta,
    BigBed,
    BigWig,
    Fasta,
    File,
    GzipBed,
    GzipVcf,
    Hic,
    PathLike,
    TwoBit,
    Vcf,
    detect,
)
from ._contact import ContactReader
from ._interval import IntervalReader
from ._sequence import SequenceReader
from ._track import TrackReader

__all__ = ["open"]


def open(
    source: File | PathLike, **options: object
) -> SequenceReader | TrackReader | IntervalReader | ContactReader:
    """The reader of ``source``'s format (a path is detected first);
    ``options`` go to that reader.

    Variants have no reader of their own: they are read over a reference
    with :class:`~omicsweave.query.ConsensusReader`.
    """
    file = source if isinstance(source, File) else detect(source)
    match file:
        case Fasta() | BgzipFasta() | TwoBit():
            return SequenceReader(file, **options)  # type: ignore[arg-type]
        case BigWig():
            return TrackReader(file, **options)  # type: ignore[arg-type]
        case BigBed() | Bed() | GzipBed():
            return IntervalReader(file, **options)  # type: ignore[arg-type]
        case Hic():
            return ContactReader(file, **options)  # type: ignore[arg-type]
        case Vcf() | GzipVcf() | Bcf():
            detail = "variants are read over a reference with ConsensusReader"
            raise InvalidError(f"{file.path}: {detail}", what=str(file.path), detail=detail)
    raise TypeError(f"no reader reads {type(file).__name__}")

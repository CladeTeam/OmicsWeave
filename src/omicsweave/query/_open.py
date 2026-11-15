"""Opening the reader a file's format calls for."""

from __future__ import annotations

from ..errors import InvalidError
from ..formats import (
    Bam,
    Bcf,
    Bed,
    BedGraph,
    BgzipFasta,
    BigBed,
    BigWig,
    Cram,
    Fasta,
    File,
    GzipBed,
    GzipBedGraph,
    GzipVcf,
    GzipWig,
    Hic,
    PathLike,
    TwoBit,
    Vcf,
    Wig,
    detect,
)
from ._alignment import AlignmentReader
from ._contact import ContactReader
from ._interval import IntervalReader, _bigbed_name
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

    VCF and BED text are not read by region: convert them to BCF and
    bigBed first.
    """
    file = source if isinstance(source, File) else detect(source)
    match file:
        case Fasta() | BgzipFasta() | TwoBit():
            return SequenceReader(file, **options)  # type: ignore[arg-type]
        case BigWig():
            return TrackReader(file, **options)  # type: ignore[arg-type]
        case BigBed():
            return IntervalReader(file, **options)  # type: ignore[arg-type]
        case Bed() | GzipBed():
            detail = (
                "BED text is not read by region; convert it to bigBed with "
                f"`ow convert {file.path} {_bigbed_name(file)} --chrom-sizes <genome.sizes>`"
            )
            raise InvalidError(f"{file.path}: {detail}", what=str(file.path), detail=detail)
        case Hic():
            return ContactReader(file, **options)  # type: ignore[arg-type]
        case BedGraph() | GzipBedGraph() | Wig() | GzipWig():
            detail = "signal text is not read by region; convert it to bigWig (ow convert)"
            raise InvalidError(f"{file.path}: {detail}", what=str(file.path), detail=detail)
        case Bam() | Cram():
            return AlignmentReader(file, **options)  # type: ignore[arg-type]
        case Bcf():
            return VariantReader(file, **options)  # type: ignore[arg-type]
        case Vcf() | GzipVcf():
            detail = "VCF text is not read by region; convert it to BCF (ow convert --index)"
            raise InvalidError(f"{file.path}: {detail}", what=str(file.path), detail=detail)
    raise TypeError(f"no reader reads {type(file).__name__}")

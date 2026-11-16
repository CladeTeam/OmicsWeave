"""Motif collections: the FFM format OmicsWeave reads, and the motif texts of
other tools that convert into it.

Motif texts are not read directly: ``ow convert`` turns one file or a list
of files (``ow convert pwms/*.txt cisbp.ffm``) into one ``.ffm``.
:class:`FfmWriter` writes one from :class:`~omicsweave.types.Motif` values;
:class:`~omicsweave.query.MotifReader` reads it.
Every text may also be gzip compressed. Values are stored as the source
gives them; count matrices (JASPAR, TRANSFAC, HOCOMOCO PCM) compress far
better than probabilities and are the recommended source.
"""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from typing import TYPE_CHECKING, Literal

from .. import _native
from ._base import File, PathLike, Writer

if TYPE_CHECKING:
    from ..types import Alphabet, Motif

__all__ = [
    "Bamm",
    "CisBp",
    "Dpcm",
    "Dpwm",
    "Ffm",
    "FfmWriter",
    "Glam2",
    "Homer",
    "Jaspar",
    "JasparPfm",
    "Meme",
    "Pcm",
    "ProBound",
    "Pwm",
    "Transfac",
]


class Ffm(File):
    """FFM (Factor Frequency Model) motif collection (``.ffm``)."""

    __slots__ = ()
    format = "ffm"


class Meme(File):
    """MEME motif text (``.meme``): ``MOTIF`` records with letter-probability
    matrices, in the alphabet the file declares."""

    __slots__ = ()
    format = "meme"


class Jaspar(File):
    """JASPAR count matrices (``.jaspar``): ``>ID name`` and ``A [ ... ]`` rows."""

    __slots__ = ()
    format = "jaspar"


class JasparPfm(File):
    """JASPAR raw PFM: ``>ID name`` and four rows of counts, A, C, G, T
    (also HOCOMOCO's JASPAR-format download)."""

    __slots__ = ()
    format = "jaspar-pfm"


class Transfac(File):
    """TRANSFAC matrix records (``AC``, ``ID``, ``P0`` lines, ``//``)."""

    __slots__ = ()
    format = "transfac"


class Homer(File):
    """HOMER motifs (``.motif``): ``>consensus name threshold`` and rows of
    probabilities."""

    __slots__ = ()
    format = "homer"


class CisBp(File):
    """CIS-BP PWM, one motif per file (``Pos A C G T`` header); the file name
    is the motif ID."""

    __slots__ = ()
    format = "cisbp"


class Pcm(File):
    """One row of counts or frequencies per position under a ``>`` header:
    HOCOMOCO ``.pcm`` / ``.pfm``, Cluster-Buster ``.cb``."""

    __slots__ = ()
    format = "pcm"


class Pwm(File):
    """One row of log weights per position (HOCOMOCO ``.pwm``)."""

    __slots__ = ()
    format = "pwm"


class Dpcm(File):
    """Dinucleotide counts, sixteen per adjacent pair (HOCOMOCO ``.dpcm``)."""

    __slots__ = ()
    format = "dpcm"


class Dpwm(File):
    """Dinucleotide log weights (HOCOMOCO ``.dpwm``)."""

    __slots__ = ()
    format = "dpwm"


class Bamm(File):
    """BaMM inhomogeneous Markov model (``.ihbcp``), one motif per file."""

    __slots__ = ()
    format = "bamm"


class ProBound(File):
    """ProBound / MotifCentral JSON model; each binding mode is a motif."""

    __slots__ = ()
    format = "probound"


class Glam2(File):
    """GLAM2 text output (``glam2.txt``); each alignment is a gapped profile
    of residue, deletion and insertion counts."""

    __slots__ = ()
    format = "glam2"


class FfmWriter(Writer[Ffm]):
    """Writes an FFM collection: motifs over ``alphabet`` in the order
    added, IDs unique, annotation values for the keys of ``schema`` (name →
    ``"string"``, ``"integer"`` or ``"float"``, in stored order; keys a
    motif lacks are stored missing). Blocks are compressed as motifs come
    and the file is written on :meth:`close`.

    ``block_records`` and ``block_values`` bound a block's motifs and table
    cells: smaller blocks make single lookups cheaper, larger ones compress
    better.
    """

    def __init__(
        self,
        path: PathLike,
        alphabet: Alphabet,
        schema: Mapping[str, Literal["string", "integer", "float"]] | None = None,
        *,
        block_records: int = 4096,
        block_values: int = 4096,
    ) -> None:
        schema = dict(schema or {})
        native = _native.FfmWriter(
            os.fspath(path), alphabet, list(schema.items()), block_records, block_values
        )
        super().__init__(Ffm(path), native)
        self.alphabet = alphabet
        self._keys = tuple(schema)

    def add(self, motif: Motif) -> None:
        """Appends ``motif``; its model is checked against the alphabet."""
        unknown = set(motif.annotation) - set(self._keys)
        if unknown:
            raise ValueError(f"annotation keys {sorted(unknown)} are not in the schema")
        values = ["" if (v := motif.annotation.get(k)) is None else str(v) for k in self._keys]
        self._native.add(motif, values)

    def add_all(self, motifs: Iterable[Motif]) -> None:
        for m in motifs:
            self.add(m)

    def __len__(self) -> int:
        return self._native.__len__()

    def _finish(self) -> tuple[File, ...]:
        self._native.close()
        return ()

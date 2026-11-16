"""Results printed as MEME Suite's tools print them, byte for byte.

Each function takes a result and what the text shows beside it: motif IDs
(and alternate IDs) by motif number, sequence names by sequence number,
and, where matched letters are printed, the scanned sequences. ``origins``
are added to every printed position (a region's start, say).
"""

from __future__ import annotations

from collections.abc import Sequence as Seq

import numpy as np

from .. import _native
from ..types import Alphabet, StringVector
from ._profile import ProfileAlignments
from ._reduce import ClusterTable, HitTable, MastReport
from ._sequences import Sequences, as_strings

__all__ = ["fimo", "glam2scan", "mast_hit_list", "mcast"]


def _names(names: StringVector | Seq[str]) -> StringVector:
    return names if isinstance(names, StringVector) else StringVector.from_list(names)


def _origins(origins: Seq[int] | np.ndarray | None) -> np.ndarray | None:
    return None if origins is None else np.ascontiguousarray(origins, dtype=np.uint64)


def fimo(
    hits: HitTable,
    *,
    alphabet: Alphabet,
    motif_ids: Seq[str],
    names: StringVector | Seq[str],
    sequences: Sequences,
    alt_ids: Seq[str] | None = None,
    origins: Seq[int] | np.ndarray | None = None,
) -> str:
    """FIMO 5.5.9 ``--text --no-qvalue``: a header, then one line per hit
    with 1-based inclusive positions and the matched letters (the
    alphabet's complement symbols on the reverse strand)."""
    n = _names(names)
    s = as_strings(sequences)
    return _native.text_fimo(
        hits,
        alphabet,
        list(motif_ids),
        list(alt_ids or [""] * len(motif_ids)),
        (n.offsets, n.data),
        (s.offsets, s.data),
        _origins(origins),
    )


def mast_hit_list(
    report: MastReport,
    *,
    motif_ids: Seq[str],
    names: StringVector | Seq[str],
    alt_ids: Seq[str] | None = None,
    numbers: Seq[int] | None = None,
    stranded: bool = True,
) -> str:
    """The lines of ``mast -hit_list`` after its two ``#`` header lines:
    ``name ±number id alt start end score p-value`` per tiling window.
    ``numbers`` are the motifs' 1-based numbers in their source file
    (their order by default); ``stranded=False`` prints no strand."""
    n = _names(names)
    numbers = list(numbers) if numbers is not None else list(range(1, len(motif_ids) + 1))
    return _native.text_mast_hit_list(
        report,
        list(motif_ids),
        list(alt_ids or [""] * len(motif_ids)),
        numbers,
        (n.offsets, n.data),
        stranded,
    )


def mcast(
    clusters: ClusterTable,
    *,
    names: StringVector | Seq[str],
    sequences: Sequences,
    origins: Seq[int] | np.ndarray | None = None,
    evalue: float | None = None,
    pvalue: float | None = None,
    qvalue: float | None = None,
) -> str:
    """MCAST's ``--text`` table. With statistics only clusters within all
    three thresholds print; unset ones take MCAST's defaults (``E <= 10``,
    no p- or q-value limit). ``--output-pthresh x`` is ``pvalue=x,
    evalue=inf``; ``--output-qthresh x`` is ``qvalue=x, evalue=1``. Every
    cluster prints without statistics."""
    n = _names(names)
    s = as_strings(sequences)
    return _native.text_mcast(
        clusters,
        (n.offsets, n.data),
        (s.offsets, s.data),
        _origins(origins),
        (evalue, pvalue, qvalue),
    )


def glam2scan(
    alignments: ProfileAlignments,
    *,
    alphabet: Alphabet,
    names: StringVector | Seq[str],
    sequences: Sequences,
    header: tuple[str, str] | None = None,
) -> str:
    """``glam2scan -t``'s alignments; ``header`` is ``(version, command
    line)`` for the three header lines before them."""
    n = _names(names)
    s = as_strings(sequences)
    return _native.text_glam2(
        alignments, alphabet, (n.offsets, n.data), (s.offsets, s.data), header
    )

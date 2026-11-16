"""MEME Suite's motif scanners as compositions of the framework: each is a
scoring preset, one reduction and one scan."""

from __future__ import annotations

import dataclasses
from collections.abc import Sequence as Seq
from typing import Any, Literal

import numpy as np

from ..types import Alphabet, StringVector
from ._background import Background
from ._defaults import default
from ._profile import ProfileAlignments, ProfileScoring, compile_profiles
from ._program import MotifSource, compile, opened
from ._reduce import Clusters, ClusterTable, Hits, HitTable, Mast, MastReport
from ._scoring import Scoring, Strands
from ._sequences import Sequences
from .stats import Adjust, BestOf, ClusterStatistics, Combine

__all__ = ["fimo", "glam2scan", "mast", "mcast"]


def _given(**options: Any) -> dict[str, Any]:
    return {k: v for k, v in options.items() if v is not None}


def fimo(
    motifs: MotifSource,
    sequences: Sequences,
    *,
    background: Background,
    pvalue: float = 1e-4,
    strands: Strands = "both",
    alphabet: Alphabet | None = None,
    threads: int | None = None,
) -> HitTable:
    """FIMO: every window with p-value at most ``pvalue`` under FIMO's
    scoring (:meth:`Scoring.fimo`). q-values are a separate step
    (``stats.BenjaminiHochberg()(hits.pvalue)``)."""
    program = compile(
        motifs, background, Scoring.fimo(), strands=strands, alphabet=alphabet, threads=threads
    )
    return program.scan(sequences, Hits(pvalue=pvalue), threads=threads)


def mast(
    motifs: MotifSource,
    sequences: Sequences,
    *,
    background: Background,
    letters: Seq[float] | np.ndarray | None = None,
    mode: Literal["combine", "separate", "forward"] = "combine",
    hit_threshold: float | None = None,
    evalue_threshold: float | None = None,
    best_of: BestOf | None = None,
    combine: Combine | None = None,
    adjust: Adjust | None = None,
    alphabet: Alphabet | None = None,
    threads: int | None = None,
) -> MastReport:
    """MAST: sequences scored by the combined p-value of their motifs' best
    windows, E-values over the sequences, tilings of strong windows.
    ``background`` gives the scores and ``letters`` the p-value letter
    frequencies (MAST's built-in ones when not given). An alphabet without
    complement is scanned forward only. Options left ``None`` take MAST's
    defaults (:class:`Mast`)."""
    with opened(motifs, alphabet) as (motifs, alphabet):
        if alphabet.complement is None:
            mode = "forward"
        given = None if letters is None else np.ascontiguousarray(letters, dtype=np.float64)
        scoring = Scoring.mast(alphabet, letters=given, mode=mode)
        strands: Strands = "+" if mode == "forward" else "both"
        program = compile(
            motifs, background, scoring, strands=strands, alphabet=alphabet, threads=threads
        )
    reduction = Mast(
        mode,
        **_given(
            hit_threshold=hit_threshold,
            evalue_threshold=evalue_threshold,
            best_of=best_of,
            combine=combine,
            adjust=adjust,
        ),
    )
    return program.scan(sequences, reduction, threads=threads)


def mcast(
    motifs: MotifSource,
    sequences: Sequences,
    *,
    background: Background,
    pthresh: float | None = None,
    max_gap: int | None = None,
    egcost: float | None = None,
    threshold: float | None = None,
    statistics: ClusterStatistics | bool = True,
    seed: int = 0,
    alphabet: Alphabet | None = None,
    threads: int | None = None,
) -> ClusterTable:
    """MCAST: clusters of motif hits along each sequence, scored as MCAST
    scores them (:meth:`Scoring.mcast`, both strands). ``statistics=True``
    computes MCAST's p-, E- and q-values with ``seed``; a
    :class:`~omicsweave.motif.stats.ClusterStatistics` replaces any of
    their parts; ``False`` reports clusters only. Options left ``None``
    take MCAST's defaults (:class:`Clusters`)."""
    reduction = Clusters(
        **_given(pthresh=pthresh, max_gap=max_gap, egcost=egcost), threshold=threshold
    )
    if statistics is True:
        statistics = ClusterStatistics.mcast(seed=seed, max_gap=reduction.max_gap)
    if statistics is not False:
        reduction = dataclasses.replace(reduction, statistics=statistics)
    with opened(motifs, alphabet) as (motifs, alphabet):
        scoring = Scoring.mcast(background, alphabet)
        program = compile(
            motifs, background, scoring, strands="both", alphabet=alphabet, threads=threads
        )
    return program.scan(sequences, reduction, threads=threads)


def glam2scan(
    profiles: MotifSource,
    sequences: Sequences,
    *,
    background: Background | None = None,
    scoring: ProfileScoring | None = None,
    top: int | None = None,
    names: StringVector | Seq[str] | None = None,
    strands: Strands = "both",
    alphabet: Alphabet | None = None,
    threads: int | None = None,
) -> ProfileAlignments:
    """GLAM2Scan: each gapped profile's best ``top`` alignments (GLAM2's 25
    by default) over ``sequences``, ties broken by ``names``.
    ``background`` and ``scoring`` default to GLAM2's for the alphabet."""
    with opened(profiles, alphabet) as (profiles, alphabet):
        background = Background.glam2(alphabet) if background is None else background
        program = compile_profiles(
            profiles, background, scoring, strands=strands, alphabet=alphabet
        )
    top = default("glam2scan", "top")() if top is None else top
    return program.scan(sequences, top=top, names=names, threads=threads)

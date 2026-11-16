"""Motif matching: FIMO, MAST, MCAST and GLAM2Scan as one composable
toolkit over FFM motif collections.

A scan is built in three steps, each explicit:

1. **Models and policy.** Motifs come from an FFM collection
   (:class:`~omicsweave.query.MotifReader`) or as
   :class:`~omicsweave.types.Motif` values. A :class:`Background` is always
   given, never inferred: from frequencies, a MEME background file, a
   preset or sequences, of any Markov order. A :class:`Scoring` holds every
   rule that turns frequencies into scores and exact p-values (presets
   :meth:`~Scoring.fimo`, :meth:`~Scoring.exact`, :meth:`~Scoring.mast`,
   :meth:`~Scoring.mcast`).
2. **Compile.** :func:`compile` turns them into a :class:`Program`; several
   programs, under different policies, :func:`combine` into one.
3. **Scan.** :meth:`Program.scan` reads each sequence once for every motif
   and keeps what its reductions ask for, at the level each reduces to:
   windows (:class:`Hits`, :class:`Scores`), sequences (:class:`Best`,
   :class:`Counts`, :class:`Clusters`) or the collection (:class:`Mast`).
   Gapped profiles compile with :func:`compile_profiles` and align with
   :meth:`ProfileProgram.scan`.

Every piece is usable on its own: :meth:`Program.model` opens a compiled
motif (:class:`CompiledMotif`: log-odds :class:`Term` tables, their integer
form, the exact :class:`NullDistribution`) for thresholds and statistics
built outside a scan; backgrounds estimate, symmetrize and sample; the
statistics steps run on any numpy array.

Sequence-level and collection-level statistics (best-of-n, combined
p-values, E- and q-values, tail fits) are replaceable steps in
:mod:`~omicsweave.motif.stats`; MEME's formulas are one choice each, and a
Python function of the right shape works wherever a step is taken.

:func:`fimo`, :func:`mast`, :func:`mcast` and :func:`glam2scan` run MEME
Suite 5.5.9's tools as compositions of the above and match their output;
:mod:`~omicsweave.motif.text` prints results as those tools do.

Sequences are given in memory as ASCII bytes (``str``, ``bytes``,
``uint8`` arrays, :class:`~omicsweave.types.Sequence` or a
:class:`~omicsweave.types.StringVector` of many); reading FASTA and naming
sequences is the caller's. Results do not depend on the thread count.

>>> from omicsweave import motif, types
>>> dna = types.Alphabet.dna()
>>> m = types.Motif.from_pfm("TATA", [[0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0]])
>>> program = motif.compile([m], motif.Background.uniform(dna), alphabet=dna)
>>> hits = program.scan("GGTATAGG", motif.Hits(pvalue=0.01))
>>> hits.start.tolist(), hits.strand.tolist()
([2, 2], [1, -1])
"""

from ..types import Alphabet, Factor, Gaps, Motif
from . import stats, text
from ._background import Background
from ._compiled import CompiledMotif, NullDistribution, Term
from ._profile import (
    DirichletMixture,
    ProfileAlignments,
    ProfileProgram,
    ProfileScoring,
    compile_profiles,
)
from ._program import MotifSource, Part, Program, combine, compile
from ._reduce import (
    Best,
    BestTable,
    Clusters,
    ClusterTable,
    CountMatrix,
    Counts,
    Hits,
    HitTable,
    Mast,
    MastHits,
    MastMotifs,
    MastReport,
    Reduction,
    Scores,
    ScoreTracks,
)
from ._scoring import (
    Ambiguity,
    AmbiguityRule,
    ColumnPrior,
    LogScale,
    Quantizer,
    Scoring,
    Strands,
    SubstitutionPrior,
)
from ._sequences import Sequences
from ._tools import fimo, glam2scan, mast, mcast

__all__ = [
    "Alphabet",
    "Ambiguity",
    "AmbiguityRule",
    "Background",
    "Best",
    "BestTable",
    "ClusterTable",
    "Clusters",
    "ColumnPrior",
    "CompiledMotif",
    "CountMatrix",
    "Counts",
    "DirichletMixture",
    "Factor",
    "Gaps",
    "HitTable",
    "Hits",
    "LogScale",
    "Mast",
    "MastHits",
    "MastMotifs",
    "MastReport",
    "Motif",
    "MotifSource",
    "NullDistribution",
    "Part",
    "ProfileAlignments",
    "ProfileProgram",
    "ProfileScoring",
    "Program",
    "Quantizer",
    "Reduction",
    "ScoreTracks",
    "Scores",
    "Scoring",
    "Sequences",
    "Strands",
    "SubstitutionPrior",
    "Term",
    "combine",
    "compile",
    "compile_profiles",
    "fimo",
    "glam2scan",
    "mast",
    "mcast",
    "stats",
    "text",
]

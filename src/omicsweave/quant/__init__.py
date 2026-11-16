"""Quantifying aligned reads: one pass over coordinate-sorted BAMs into
coverage tracks, count matrices, single-cell ATAC fragments, methylation,
Hi-C contacts and allele counts, in parallel
within a budget of threads and memory.

- :class:`Coverage`: bigWig of aligned bases, spans, read ends or
  fragments, both strands or split, raw or scaled.
- :class:`Counts` (presets :meth:`~Counts.rna_single_cell`,
  :meth:`~Counts.rna_bulk`, :meth:`~Counts.atac`, :meth:`~Counts.chip`):
  rows (barcodes, read groups, input files) × features (:class:`GeneExons`,
  :class:`GeneSpans`, :class:`Promoters`, :class:`GeneActivity`,
  :class:`Peaks`, :class:`Bins`, or any interval set as :class:`Regions`)
  of reads, molecules, read pairs, fragments or Tn5 cut sites, as a 10x
  directory with a CellPack (``.cpk``), Matrix Market or CSV matrix.
- :class:`Atac`: ``fragments.tsv.gz`` with its tabix index, per-barcode
  metrics, and peak / tile / gene-activity matrices.
- :class:`Junctions`: splice junctions with STAR's ``SJ.out.tab`` columns,
  and a junction × barcode (or sample) matrix (:class:`JunctionMatrix`).
- :class:`Methylation`: per-cytosine table, beta / coverage bigWigs, M-bias
  and region × cell matrices (:class:`MethylMatrix`) from Bismark ``XM``,
  ``MM`` / ``ML`` or bisulfite reads against the reference.
- :class:`Contacts`: Hi-C 4DN pairs and ``.hic`` with per-cell statistics.
- :class:`Alleles`: variant × barcode (or sample) ``ref`` / ``alt`` /
  ``other`` matrices at SNVs, per-variant table, gene-level
  allele-specific expression (:class:`AlleleGenes`).
- Gene matrices with ``velocity`` (spliced / unspliced / ambiguous layers)
  and ``labeled`` (:class:`Labeled`) layers; :class:`Labeling`: per-read
  T>C conversion histograms and mismatch rates (SLAM-seq, 4sU).

A :class:`Task` names the BAMs, the products and the output directory;
:func:`plan` checks it against the budget and returns the :class:`Plan`
whose :meth:`Plan.run` makes the pass; :func:`run` does both. Defaults
and presets come from the native layer, which also checks every option. BAMs
must be sorted by coordinate (``ow convert in.bam sorted.bam --sort``).
The native layer does the work; the results do not depend on the number of
threads.

Molecule deduplication encodes uppercase ACGT UMIs up to 31 bases without
collisions. Other UMIs use 63-bit fingerprints; read-name-based molecule
keys use 64-bit fingerprints. Fingerprint collisions are not disambiguated.
Matching the two mates of a pair compares their full read names separately.
"""

from .._budget import Budget, MemoryPlan
from ._features import Bins, Features, GeneActivity, GeneExons, GeneSpans, Peaks, Promoters, Regions
from ._products import (
    AlleleGenes,
    Alleles,
    Atac,
    Contacts,
    Context,
    Counts,
    Coverage,
    Filter,
    JunctionMatrix,
    Junctions,
    Labeled,
    Labeling,
    Methylation,
    MethylMatrix,
    Product,
    Strandedness,
)
from ._run import Input, OutputFile, Plan, Report, Task, plan, run

__all__ = [
    "AlleleGenes",
    "Alleles",
    "Atac",
    "Bins",
    "Budget",
    "Contacts",
    "Context",
    "Counts",
    "Coverage",
    "Features",
    "Filter",
    "GeneActivity",
    "GeneExons",
    "GeneSpans",
    "Input",
    "JunctionMatrix",
    "Junctions",
    "Labeled",
    "Labeling",
    "MemoryPlan",
    "MethylMatrix",
    "Methylation",
    "OutputFile",
    "Peaks",
    "Plan",
    "Product",
    "Promoters",
    "Regions",
    "Report",
    "Strandedness",
    "Task",
    "plan",
    "run",
]

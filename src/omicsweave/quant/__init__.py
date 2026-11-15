"""Quantifying aligned reads: one pass over coordinate-sorted BAMs into
coverage tracks, count matrices, single-cell ATAC fragments, methylation,
Hi-C contacts and allele counts, in parallel
within a budget of threads and memory.

- :class:`Coverage`: bigWig of aligned bases, spans, read ends or
  fragments, both strands or split, raw or scaled.
- :class:`Counts` (and the :func:`Genes` presets): rows (barcodes, read
  groups, input files) × features (genes, gene windows, bigBed peaks, genome
  bins) of reads, molecules, read pairs, fragments or Tn5 cut sites, as a
  10x directory with a CellPack (``.cpk``), Matrix Market or CSV matrix.
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
:func:`plan` checks it against the budget, :func:`run` makes the pass. BAMs
must be sorted by coordinate (``ow convert in.bam sorted.bam --sort``).
The native layer does the work; the results do not depend on the number of
threads.
"""

from ._products import (
    AlleleGenes,
    Alleles,
    Atac,
    Contacts,
    Counts,
    Coverage,
    Filter,
    Genes,
    JunctionMatrix,
    Junctions,
    Labeled,
    Labeling,
    Methylation,
    MethylMatrix,
    Product,
)
from ._run import MemoryPlan, Plan, Report, Task, plan, run

__all__ = [
    "AlleleGenes",
    "Alleles",
    "Atac",
    "Contacts",
    "Counts",
    "Coverage",
    "Filter",
    "Genes",
    "JunctionMatrix",
    "Junctions",
    "Labeled",
    "Labeling",
    "MemoryPlan",
    "MethylMatrix",
    "Methylation",
    "Plan",
    "Product",
    "Report",
    "Task",
    "plan",
    "run",
]

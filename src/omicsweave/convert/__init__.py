"""Converting files with a budget of threads and memory.

- FASTA (any storage) to bgzip FASTA, with its ``.fai`` and ``.gzi``
- FASTA (any storage) to ``.2bit``
- VCF (any storage) to ``.bcf``, with a ``.csi`` given ``index``
- SAM (any storage) to ``.bam``, with a ``.bai`` / ``.csi`` given ``index``
- ``.bam`` to a coordinate-sorted ``.bam`` given ``sort``, indexed given ``index``
- ``.cram`` to ``.bam`` (given the reference it was compressed against)
- BED or bedGraph (any storage) to ``.bigbed`` / ``.bb``
- bedGraph or wig (any storage) to ``.bigwig`` / ``.bw``
- ``.cool`` / ``.mcool`` to ``.hic``
- 4DN pairs, Juicer short or HiC-Pro valid pairs (any storage) to ``.hic``
- bedMethyl (any storage) to ``.bigbed`` (every column) or ``.bigwig`` (one
  column)
- MAF (any storage) to bigMaf (a ``.bigbed`` / ``.bb`` / ``.bigmaf``: one
  record per alignment block, as UCSC's ``mafToBigMaf`` and ``bedToBigBed``
  write it)
- Motif texts or FFM, from one or more files of the same format, to one
  ``.ffm`` collection; these inputs and models are loaded in memory

"Any storage" is plain text, gzip or BGZF; bigBed, bigWig and hic targets
also get a ``.sidecar`` given ``sidecar``. :func:`plan` configures a
conversion and checks it against the budget; :meth:`Plan.run` converts;
:func:`run` does both. :func:`routes` lists the pairs of formats. The
memory plan bounds streaming buffers; motif size estimates are not an RSS
cap. The native layer does the work: it tells formats apart, checks options, sizes
the buffers and writes the files, which appear only once the whole
conversion succeeded. The output does not depend on the number of
threads, nor on whether a sort had to spill to disk.
"""

from .._budget import Budget, MemoryPlan
from ._routes import routes
from ._run import ConvertWarning, Plan, Report, plan, run

__all__ = ["Budget", "ConvertWarning", "MemoryPlan", "Plan", "Report", "plan", "routes", "run"]

"""Converting a file into another format: streaming, in parallel, within a
budget of threads and memory.

- FASTA (any storage) to bgzip FASTA, with its ``.fai`` and ``.gzi``
- FASTA (any storage) to ``.2bit``
- VCF (any storage) to ``.bcf``, with a ``.csi`` given ``index``
- SAM (any storage) to ``.bam``, with a ``.bai`` / ``.csi`` given ``index``
- ``.bam`` to a coordinate-sorted ``.bam`` given ``sort``, indexed given ``index``
- BED or bedGraph (any storage) to ``.bigbed`` / ``.bb``, with its ``.sidecar``
- bedGraph or wig (any storage) to ``.bigwig`` / ``.bw``, with its ``.sidecar``
- ``.cool`` / ``.mcool`` to ``.hic``, with its ``.sidecar``
- 4DN pairs or Juicer short text (any storage) to ``.hic``, with its
  ``.sidecar``
- bedMethyl (any storage) to ``.bigbed`` (every column) or ``.bigwig`` (one
  column), with its ``.sidecar``

"Any storage" is plain text, gzip or BGZF. :func:`run` converts,
:func:`plan` only checks the budget, :func:`routes` lists the pairs of
formats. The native layer does the work: it tells formats apart, checks
options, sizes the buffers and writes the files, which appear only once
the whole conversion succeeded; this module hands it paths or file objects
and returns file objects. The output does not depend on the number of
threads, nor on whether a sort had to spill to disk.
"""

from ._plan import Budget, MemoryPlan, Plan, plan
from ._routes import routes
from ._run import ConvertWarning, Report, run

__all__ = ["Budget", "ConvertWarning", "MemoryPlan", "Plan", "Report", "plan", "routes", "run"]

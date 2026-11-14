"""Load multi-omics data for model training.

The package has three layers:

- :mod:`omicsweave.types`: the values readers take and return, whatever
  file they come from (:class:`Region`, ``Sequence``, ``ContactWindow``, ...)
- :mod:`omicsweave.formats`: one class per file format, naming a file of it
  (``BgzipFasta``, ``Fai``, ``BigWig``, ...), and the writers of each
- :mod:`omicsweave.query`: readers that look entries up by key
  (``SequenceReader``, ``ConsensusReader``, ...); :func:`open` picks one

Indexes are never built on open: ``Fasta.build_index()`` builds them, and
``FastaWriter`` saves them with what it writes. Sidecars are for loaders
only and are built explicitly (``build_sidecar()``, or ``sidecar=True`` on a
writer). :mod:`omicsweave.convert` converts a file into another format in
parallel within a memory budget (``convert.run(source, target)``); the
``ow`` command (:mod:`omicsweave.cli`) does the same from a shell.
Failures raise the exceptions of :mod:`omicsweave.errors`.

Python assembles: it names files, picks readers and shapes results. The
native layer does the work: formats, decoding, encoding, budgets.
"""

from . import convert, errors, formats, query, types
from ._native import __version__
from .errors import OmicsWeaveError
from .query import open
from .types import Region

__all__ = [
    "OmicsWeaveError",
    "Region",
    "__version__",
    "convert",
    "errors",
    "formats",
    "open",
    "query",
    "types",
]

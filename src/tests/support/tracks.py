import numpy as np

from omicsweave.formats import BigWigWriter

CHROMS = {"chr1": 50_000, "chr2": 23_000}


def track(path, rng, **kw):
    starts = np.arange(0, 40_000, 7, dtype=np.uint32)
    with BigWigWriter(path, CHROMS, **kw) as w:
        w.add("chr1", starts, starts + 5, rng.normal(size=len(starts)).astype(np.float32))
    return path

# Examples

Training-data loaders for the two starting use cases (`docs/design.md` §3),
built on a plain PyTorch `DataLoader`. They are illustrative: paths and
manifests are placeholders, and the data are taken to be converted already
(`ow convert`: bgzip FASTA, bigWig, bigBed, `.hic`; expression as a
CellPack `.cpk`).

| Example | Samples | Reads |
|---|---|---|
| [`nonb_dataloader.py`](nonb_dataloader.py) | tracks: 131 kb windows, binned signal of every track; peaks: 5 kb one-hot windows, label and per-track bin coverage | `SequenceReader`, `TrackReader`, `IntervalReader` |
| [`ceno_dataloader.py`](ceno_dataloader.py) | one 131 kb window of one cell: bp assays with valid masks and input / target split, 128 bp histone + CTCF, 2048 bp Hi-C tile, the cell's expression | `SequenceReader`, `TrackReader`, `ContactReader`, `pycellpack.Store` |

Both follow the same pattern:

- The dataset holds `omicsweave.formats` file objects (cheap, picklable) and
  opens its readers in each worker on the first sample; every later query
  reuses them.
- The sampler (or batch sampler) runs in the main process and hands out
  genomic regions, so a new epoch's windows (`set_epoch`) reach persistent
  workers.
- `prefetch_factor`, `persistent_workers` and `pin_memory` do the
  prefetching; there is no separate loader.
- CENO's batch sampler groups the cells of one window; `__getitems__` reads
  the window's bases once and the batch's expression rows in one CellPack
  decode.

Requirements beyond `omicsweave`: `torch`, and `pycellpack` (CENO).

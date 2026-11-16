# Examples

Training-data loaders for the two starting use cases ([query guide](../docs/query.md)),
built on a plain PyTorch `DataLoader`. They are illustrative: paths and
manifests are placeholders. Every file is read by region in a format built
for it (bgzip FASTA, bigWig, bigBed, `.hic`; expression as a CellPack
`.cpk`); other formats are converted once beforehand with `ow convert`
(bedGraph / wig → bigWig, BED → bigBed, `.mcool` → `.hic`).

| Example | Samples | Reads |
|---|---|---|
| [`nonb_dataloader.py`](nonb_dataloader.py) | tracks: 131 kb windows, binned signal of every track; peaks: 5 kb one-hot windows, label and per-track bin coverage | `SequenceReader`, `TrackReader`, `IntervalReader` |
| [`ceno_dataloader.py`](ceno_dataloader.py) | one 131 kb window of one cell: bp assays with valid masks and input / target split, 128 bp histone + CTCF, 2048 bp Hi-C tile, the cell's expression | `SequenceReader`, `TrackReader`, `ContactReader`, `pycellpack.Store` |

Both follow the same pattern:

- Queries return arrays directly: `fetch(...).bases`, `signal(...).values`,
  `bins(...).values` and `matrix(...).values`. Readers allocate outputs and
  manage decoding workspaces; the dataset assembles these results into its
  sample arrays. Optional `out=` reuse is described in the main README.
- The dataset holds `omicsweave.formats` file objects (cheap, picklable) and
  opens its readers in each worker on the first sample; every later query
  reuses them.
- CENO and the NonB track pipeline make one `omicsweave.cache.BlockCache`
  per worker (it is process-local and never pickled) and pass it to its
  readers with `cache=`, so windows that overlap or recur reuse the decoded
  blocks they share.
- The sampler (or batch sampler) runs in the main process and hands out
  genomic regions, so a new epoch's windows (`set_epoch`) reach persistent
  workers.
- `prefetch_factor`, `persistent_workers` and `pin_memory` do the
  prefetching; there is no separate loader.
- CENO's batch sampler groups the cells of one window; `__getitems__` reads
  the window's bases once and the batch's expression rows in one CellPack
  decode.

Requirements beyond `omicsweave`: `torch`, and `pycellpack` (CENO).

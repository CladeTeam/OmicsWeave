<h1 align="center">OmicsWeave</h1>

<p align="center">
  <strong>Genomic files, read by region, straight into model training.</strong>
</p>

OmicsWeave reads the files genomics already uses (FASTA, 2bit, bigWig,
bigBed, BED, `.hic`, VCF/BCF) by genomic region and returns NumPy arrays,
so training samples are built from the original files instead of an
offline copy. It also writes those formats and converts between them in
parallel within a memory budget. The engine is Rust; the Python API is a
thin, typed layer over it.

- **One query, one array.** Ask for a region, get its bases, signal,
  features or contact tile as NumPy columns. Fixed-size results can be
  written into your own buffers.
- **Matches the reference tools.** Every format is checked against the
  tool it replaces (htslib, samtools, bgzip, bcftools, libBigWig, lib2bit,
  straw); indexes, encodings and consensus sequences are byte-identical.
- **Built for data loaders.** Readers open once and answer any number of
  queries, release the GIL while decoding, and pickle into DataLoader
  workers.
- **Conversion at scale.** FASTA, VCF, SAM, BED and cooler files convert to
  bgzip FASTA, 2bit, BCF, BAM, bigBed and `.hic` with indexes, streaming
  and in parallel; the output does not depend on the thread count.

## Install

OmicsWeave builds from source and needs Python 3.12+, a Rust toolchain
(1.85+), a C compiler and CMake 3.26+ (HDF5 is built and linked
statically).

```sh
pip install git+https://github.com/CladeTeam/OmicsWeave
# or, in a clone
uv sync
```

The install provides the `omicsweave` package and the `ow` command.

## Quick start

```python
import omicsweave as ow
from omicsweave.formats import BgzipFasta
from omicsweave.query import SequenceReader, TrackReader

BgzipFasta("hg38.fa.gz").build_index()  # once: writes .fai and .gzi

with SequenceReader("hg38.fa.gz") as genome, TrackReader("atac.bw") as atac:
    region = ow.Region("chr1", 1_000_000, 1_131_072)  # 0-based, half-open
    bases = genome.fetch(region).bases  # uint8 ASCII, [131072]
    signal = atac.signal(region).values  # float32, [131072]; NaN = no data
```

A region may also be given as a tuple: `genome.fetch(("chr1", 0, 100))`.

## Formats

| Format | Read | Write | Convert from | Convert to |
|---|:---:|:---:|:---:|:---:|
| FASTA: plain, bgzip (`.fai`, `.gzi`) | ✓ | ✓ | ✓ (also plain gzip) | bgzip |
| 2bit | ✓ | | | ✓ |
| bigWig | ✓ | ✓ | | |
| bigBed | ✓ | ✓ | | ✓ |
| BED, BED.gz | ✓ | | ✓ | |
| `.hic` v8 / v9 | ✓ | ✓ | | ✓ |
| `.cool`, `.mcool` | | | ✓ | |
| VCF / BCF | consensus sequences | | VCF | BCF (+ `.csi`) |
| SAM / BAM | | | SAM | BAM (+ `.bai` / `.csi`) |

Formats are told apart by content, not by name: `omicsweave.formats.detect`
names the format of any of these files.

## Reading

The package has three layers: `omicsweave.formats` names files,
`omicsweave.query` reads them, and `omicsweave.types` holds the keys and
results. Every reader takes a file object or a path, works as a context
manager, and returns dataclasses of NumPy arrays that carry their key.

| Reader | Files | Key | Results |
|---|---|---|---|
| `SequenceReader` | FASTA, bgzip FASTA, 2bit | `Region` | `Sequence` |
| `TrackReader` | bigWig | `Region` | `Signal`, `SignalIntervals` |
| `IntervalReader` | bigBed, BED, BED.gz | `Region` | `Intervals` |
| `ContactReader` | `.hic` | `ContactWindow` | `ContactMatrix`, `ContactPixels` |
| `ConsensusReader` | VCF/BCF + a reference | `Region` + alleles | `ConsensusSequence`, `ConsensusBatch` |

### Sequences

```python
import numpy as np
from omicsweave.formats import BgzipFasta, Fai, Gzi

genome = SequenceReader(BgzipFasta("hg38.fa.gz"))  # indexes beside the file
genome = SequenceReader("hg38.fa.gz", Fai("idx/hg38.fai"), Gzi("idx/hg38.gzi"))
genome = SequenceReader("hg38.2bit", soft_mask=True)  # 2bit: no index needed

buffer = np.empty(131_072, np.uint8)
genome.fetch(("chr1", 0, 131_072), out=buffer)  # decode in place
genome.chroms  # {'chr1': 248956422, ...}
```

Indexes are never built on open: a missing `.fai` or `.gzi` raises
`IoError`. Build them with `build_index()`, or write the file with
`FastaWriter`, which saves them.

### Signal and features

```python
from omicsweave.query import IntervalReader

with TrackReader("H3K27ac.bw") as track:
    track.signal(region).values  # one float32 per base
    runs = track.intervals(region)  # stored runs, clipped: start, end, value columns
    track.summary  # bases covered, min, max, sum, sum of squares

with IntervalReader("peaks.bb") as peaks:  # bigBed, or BED loaded whole
    hits = peaks.fetch(region)  # overlapping features, by start
    hits.start, hits.end  # uint32 columns
    hits.annotation.to_list()  # fields after the third, tab-joined
    hits.annotation.to_arrow()  # zero-copy pyarrow LargeStringArray
```

### Contact maps

```python
from omicsweave.query import ContactReader
from omicsweave.types import ContactWindow

with ContactReader("k562.hic") as hic:
    window = ContactWindow(("chr1", 0, 5_000_000), ("chr1", 0, 5_000_000), 5_000)
    hic.matrix(window).values  # float32 [1000, 1000], raw counts
    pixels = hic.pixels(window)  # nonzero pixels: row, col, value
```

Counts are returned as stored; normalization is left to the caller.

### Personal sequences

`ConsensusReader` applies the records of a VCF or BCF to a reference, as
`bcftools consensus` 1.23 does, entirely in memory.

```python
from omicsweave.query import ConsensusOptions, ConsensusReader
from omicsweave.types import Haplotype, Samples

with ConsensusReader("cohort.vcf.gz", genome, ConsensusOptions(include="QUAL>20")) as cohort:
    cohort.fetch(region, Haplotype("NA12878", "1")).text
    cohort.fetch(region, Samples(["NA12878", "NA12891"]))  # IUPAC codes
    batch = cohort.fetch_batch(region, [Haplotype(s, h) for s in cohort.samples for h in "12"])
```

A batch reads the reference and finds the records of its region once, and
returns every sequence packed in one array.

## Writing

```python
from omicsweave.formats import BigBedSchema, BigBedWriter, BigWigWriter, FastaWriter

with FastaWriter("out.fa.gz") as fasta:  # bgzip; saves .fai and .gzi
    fasta.add("chr1", "ACGT...")

with BigWigWriter("out.bw", {"chr1": 248_956_422}) as bw:
    bw.add("chr1", starts, ends, values)  # columns or scalars

with BigBedWriter("out.bb", chroms, BigBedSchema(field_count=6)) as bb:
    bb.add("chr1", starts, ends, ["name\t0\t+", ...])
```

`HicWriter` writes binned maps and `HicBuilder` bins contact pairs at
every resolution, as Juicer's `pre` does. Writers finish on leaving the
`with` block and return the companion files they wrote.

## Converting

```sh
ow convert hg38.fa.gz hg38.2bit
ow convert cohort.vcf.gz cohort.bcf --index -t 16 -m 8G
ow convert reads.sam reads.bam --index
ow convert peaks.bed peaks.bb --chrom-sizes hg38.fa.fai
ow convert k562.mcool k562.hic --resolutions 1000,5000
ow convert genome.txt genome.out --from fasta --to bgzip-fasta
```

The source's format is read from its content and the target's from its
name; `--from` and `--to` state them. `--dry-run` shows the memory plan
without converting, and `ow convert --help` lists every conversion. The
target and its indexes appear only once the whole conversion succeeded.

The same from Python:

```python
report = ow.convert.run("cohort.vcf.gz", "cohort.bcf", index=True, threads=16, memory="8G")
report.target, report.companions  # Bcf('cohort.bcf'), (Csi('cohort.bcf.csi'),)
ow.convert.plan("hg38.fa.gz", "hg38.2bit", memory="1G").memory
```

Sorts that outgrow the memory budget spill to disk beside the target (or
to `tmp_dir`); the output is the same either way.

## Training with PyTorch

Readers fit a plain `DataLoader`: keep the file objects in the dataset,
open readers in the worker on first use, and let the sampler hand out
regions.

```python
from omicsweave.formats import BigWig
from torch.utils.data import DataLoader, Dataset


class Windows(Dataset):
    def __init__(self, reference, track):
        self.reference, self.track, self.readers = reference, track, None

    def __getitem__(self, region):
        if self.readers is None:  # once per worker
            self.readers = SequenceReader(self.reference), TrackReader(self.track)
        genome, signal = self.readers
        return genome.fetch(region).bases, signal.signal(region).values


loader = DataLoader(
    Windows(BgzipFasta("hg38.fa.gz"), BigWig("atac.bw")),
    sampler=windows,
    batch_size=8,
    num_workers=16,
    prefetch_factor=4,
    persistent_workers=True,
)
```

[`examples/`](examples/) has complete loaders for two use cases: NonB
(binned multi-track signal; peak classification) and CENOMultiome (per-cell
bp and coarse tracks, Hi-C tiles and expression read with CellPack).

## Behavior

- Positions are 0-based and half-open and fit in `uint32`.
- Failures raise the exceptions of `omicsweave.errors`, each also derived
  from the closest built-in (`UnknownChromError` is a `LookupError`,
  `IoError` an `OSError`, ...).
- Readers release the GIL while decoding, so threads may share one; they
  pickle as the files and options they were opened with.
- Sidecars (`.sidecar`) describe a file's data blocks for loaders; reading
  by region does not need them. Build one with `build_sidecar()` or
  `sidecar=True` on a writer.

## Fidelity

Every format is tested against the tool it replaces, built from source in
`refs/` and called in process: `samtools faidx` / `bgzip` indexes, htslib
BCF and BAM encoding and their `.csi` / `.bai` indexes, `bcftools
consensus` 1.23, libBigWig, lib2bit and straw. [`docs/design.md`](docs/design.md)
records the design and the few deliberate differences.

## Building from source

```sh
uv sync --locked                                     # build and install, editable
uv run maturin develop --release --skip-install      # rebuild after Rust changes

cargo test --workspace
uv run pytest
```

| Path | Contents |
|---|---|
| `crates/` | Rust crates `omicsweave-<name>`: one per format family, `convert`, and the `python` extension |
| `src/omicsweave/` | the Python package (`types`, `formats`, `query`, `convert`, `cli`) |
| `refs/` | reference implementations for differential tests |
| `examples/` | DataLoader examples |
| `docs/` | design notes |

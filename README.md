<h1 align="center">OmicsWeave</h1>

<p align="center">
  <strong>Genomic files, read by region, straight into model training.</strong>
</p>

OmicsWeave reads the files genomics already uses (FASTA, 2bit, bigWig,
bigBed, `.hic`, BCF, BAM, CRAM) by genomic region and returns NumPy arrays,
so training samples are built from the original files instead of an
offline copy. It also writes those formats and converts between them in
parallel within a memory budget. The engine is Rust; the Python API is a
thin, typed layer over it.

- **One query, one array.** Ask for a region, get its bases, signal,
  features, reads or contact tile as NumPy columns. Readers allocate the
  result arrays for you.
- **Matches the reference tools.** Every format is checked against the
  tool it replaces (htslib, samtools, bgzip, bcftools, libBigWig, lib2bit,
  straw); indexes, encodings and consensus sequences are byte-identical.
- **Built for data loaders.** Readers open once and answer any number of
  queries, release the GIL while decoding, and pickle into DataLoader
  workers.
- **Conversion at scale.** FASTA, VCF, SAM, CRAM, BED, bedGraph, wig, bedMethyl,
  cooler and pairs files convert to bgzip FASTA, 2bit, BCF, BAM, bigBed, bigWig and `.hic`
  with indexes, streaming and in parallel; the output does not depend on
  the thread count.

See the [documentation](docs/README.md) for installation, query, conversion,
annotation and quantification guides, and [scripts/](scripts/README.md) for benchmarks.

## Install

OmicsWeave builds from source and needs Python 3.12+, a Rust toolchain
(1.89+), a C compiler and CMake 3.26+ (HDF5 is built and linked
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
| bigWig | ✓ | ✓ | | ✓ |
| bigBed | ✓ | ✓ | | ✓ |
| BED, BED.gz | | | ✓ | |
| bedGraph, bedGraph.gz | | | ✓ | |
| wig, wig.gz | | | ✓ | |
| bedMethyl (modkit, ENCODE), bedMethyl.gz | | | ✓ | |
| `.hic` v8 / v9 | ✓ | ✓ | | ✓ |
| `.cool`, `.mcool` | | | ✓ | |
| 4DN pairs, Juicer short text (plain, gzip, bgzip) | | | ✓ | |
| VCF / BCF (`.csi`) | BCF by region (genotypes, consensus sequences) | | VCF | BCF (+ `.csi`) |
| SAM / BAM (`.bai`, `.csi`) | BAM by region | | SAM; BAM (sort) | BAM (+ `.bai` / `.csi`) |
| CRAM 3.0 / 3.1 (`.crai`) | by region (+ a reference) | | to BAM | |

Formats are told apart by content, not by name: `omicsweave.formats.detect`
names the format of any of these files.

## Reading

The package has three layers: `omicsweave.formats` names files,
`omicsweave.query` reads them, and `omicsweave.types` holds the keys and
results. Every reader takes a file object or a path, works as a context
manager, and returns dataclasses of NumPy arrays that carry their key.
Queries allocate their results automatically; no buffer sizing or separate
read/write step is needed. Results remain valid after later queries or after
the reader closes. Block caching is optional and requires no preparation
before an ordinary query.

| Reader | Files | Key | Results |
|---|---|---|---|
| `SequenceReader` | FASTA, bgzip FASTA, 2bit | `Region` / `Spliced` | `Sequence` |
| `TrackReader` | bigWig | `Region` / `Spliced` | `Signal`, `SignalBins`, `SignalIntervals` |
| `StrandedTrackReader` | a plus and a minus bigWig | stranded `Region` / `Spliced` | `Signal`, `SignalBins`, `SignalIntervals` |
| `IntervalReader` | bigBed | `Region` | `Intervals` |
| `ContactReader` | `.hic` | `ContactWindow` | `ContactMatrix`, `ContactPixels` |
| `ConsensusReader` | BCF + `.csi` + a reference | stranded `Region` + alleles | `ConsensusSequence`, `ConsensusBatch` |
| `AlignmentReader` | BAM + `.bai` / `.csi`; CRAM + `.crai` | `Region` | `Alignments` |
| `VariantReader` | BCF + `.csi` | `Region` + samples | `Genotypes` |

### Sequences

```python
from omicsweave.formats import BgzipFasta, Fai, Gzi

genome = SequenceReader(BgzipFasta("hg38.fa.gz"))  # indexes beside the file
genome = SequenceReader("hg38.fa.gz", Fai("idx/hg38.fai"), Gzi("idx/hg38.gzi"))
genome = SequenceReader("hg38.2bit", soft_mask=True)  # 2bit: no index needed

bases = genome.fetch(("chr1", 0, 131_072)).bases  # uint8 ASCII
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
    bins = track.bins(region, 128, ("mean", "coverage"))  # [2, 1024], one pass
    bins["mean"], bins["coverage"]  # mean of covered bases; fraction covered
    track.bins(region, 128, "max", phase=region.start % 128)  # bins on the genome grid
    track.bins(region, 10_000, "mean", exact=False)  # from a zoom level: an estimate
    runs = track.intervals(region)  # stored runs, clipped: start, end, value columns
    track.summary  # bases covered, min, max, sum, sum of squares

with IntervalReader("peaks.bb") as peaks:  # BED: ow convert peaks.bed peaks.bb --chrom-sizes ...
    hits = peaks.fetch(region)  # overlapping features, by start
    hits.start, hits.end  # uint32 columns
    hits.annotation.to_list()  # fields after the third, tab-joined
    hits.annotation.to_arrow()  # zero-copy pyarrow LargeStringArray
```

Bins are summarized in the reader from the stored data, without a
per-base array: `mean`, `min`, `max`, `sum`, `std` and `coverage`, as
libBigWig's `bwStats` computes them (bit for bit on full bins). Empty bins
read NaN, except `coverage`, which reads 0; the mean with missing bases as 0
is `sum / resolution`. `phase` shifts the bin grid, so the first and last
bins may be partial.

### Strands and transcripts

A `Region` has a strand (`"+"`, `"-"`, or `"."` by default), and a
`Spliced` key is several spans of a chromosome read as one sequence. Dense
results (bases, per-base signal, bins, contact tiles) run 5' to 3' of the
key's strand: on `"-"` they are reversed, and bases complemented. Sparse
results (features, runs, reads, records, pixels) keep genomic coordinates.

```python
from omicsweave.query import StrandedTrackReader
from omicsweave.types import Region, Spliced

genome.fetch(("chr1", 100, 200, "-")).text  # reverse complement
track.signal(Region("chr1", 100, 200).with_strand("-")).values  # reversed

mrna = ann.locus("ENST00000269305")  # Spliced: exons, on the transcript's strand
genome.fetch(mrna).text  # the spliced transcript, 5' to 3'
ann.locus("ENST00000269305", "cds")  # also "utr5", "utr3"; cds includes the stop codon
keys = ann.loci(coding.transcripts, "exons")  # many transcripts in one call
track.bins(mrna, 32, "mean")  # signal in transcript coordinates
mrna.to_genome([0, 1]), mrna.to_local([7_687_490])  # position <-> genome

rna = StrandedTrackReader("rna.plus.bw", "rna.minus.bw")
rna.signal(mrna).values  # the transcript's own strand
```

Every block a key's spans overlap is read and decoded once.

### Contact maps

```python
from omicsweave.query import ContactReader
from omicsweave.types import ContactWindow

with ContactReader("k562.hic") as hic:
    window = ContactWindow(("chr1", 0, 5_000_000), ("chr1", 0, 5_000_000), 5_000)
    hic.matrix(window).values  # float32 [1000, 1000], raw counts
    pixels = hic.pixels(window)  # nonzero pixels: row, col, value
```

Counts are returned as stored; normalization is left to the caller. A
region on `"-"` reverses its axis of the matrix (reverse-complement
augmentation flips both).

### Optional output reuse

For a training loop that already owns a batch array, fixed-size queries also
accept `out=`. This skips allocating a result array and writes directly into
the supplied storage; the returned result references that same array.

```python
import numpy as np

batch = np.empty((len(tracks), len(region)), dtype=np.float32)
for row, track in zip(batch, tracks, strict=True):
    result = track.signal(region, out=row)
    assert result.values is row
```

The output must be writable and C-contiguous, with the required shape and
dtype: `uint8` for bases, `float32` for signal, bins and contact matrices.
Reusing or modifying this storage also changes any result referencing it.
Readers manage decoding workspaces internally in both modes.

### Alignments

```python
from omicsweave.formats import Bam
from omicsweave.query import AlignmentReader

Bam("reads.bam").build_index()  # once, if it has none: .bai (or .csi), as samtools index

with AlignmentReader("reads.bam") as reads:  # finds .csi / .bai as htslib does
    hits = reads.fetch(region)  # start, end, flag, mapq, mate_tid, mate_start, tlen
    hits = reads.fetch(region, ("name", "seq", "cigar"))  # or "all"; also "qual", "aux"
    hits.seqs.to_list(), hits.cigars.string(0)
```

A region returns exactly the records htslib's region iterator returns, in
file order, as stored: nothing is filtered or counted.

CRAM reads the same way, given the reference it was compressed against:

```python
from omicsweave.formats import Cram

Cram("reads.cram").build_index()  # once: .crai, as samtools index

genome = SequenceReader("hg38.fa.gz")
with AlignmentReader("reads.cram", reference=genome) as reads:  # finds .crai beside it
    hits = reads.fetch(region, "all")  # what the BAM htslib writes from it holds
```

Records come as htslib decodes them, MD and NM generated where it
generates them, and each slice's reference is checked against its MD5
(`verify_checksums=False`, `decode_md=False` and `name_prefix=` change
that, as htslib's options do).
Only the blocks a query needs are decompressed: fixed fields, CIGARs,
names and qualities read without a reference (`AlignmentReader("reads.cram")`);
CRAMs with an embedded reference or none need no reference at all.

### Genotypes

```python
from omicsweave.formats import Bcf
from omicsweave.query import VariantReader

ow.convert.run("cohort.vcf.gz", "cohort.bcf", index=True)  # once: BCF + .csi
Bcf("other.bcf").build_index()  # a sorted BCF from elsewhere: .csi, as bcftools index

with VariantReader("cohort.bcf") as cohort:  # finds cohort.bcf.csi
    g = cohort.genotypes(region, ["NA12878", "NA12891"])  # default: every sample
    g.start, g.end, g.ref.to_list(), g.alt.to_list()  # one entry per record
    g.gt  # int8 [records, samples, ploidy]: allele index, -1 missing, -2 absent
    g.phased  # bool [records, samples]

    v = cohort.fetch(region, info=["DP", "AF"], format=["AD"], include="QUAL>30")
    v.genotypes.gt, v.qual, v.filters.to_list()  # the records the expression keeps
    v.info["AF"]  # float32 [records, width]; v.format["AD"]: int32 [records, samples, width]
```

The records of a region are those htslib's `bcf_itr_queryi` returns.
`genotypes` decodes only the fixed fields and `GT`; `fetch` also decodes
QUAL, FILTER and the INFO / FORMAT fields asked for, in their declared
types, padded to the widest record with htslib's missing and vector-end
values, and filters with bcftools expressions (`include` / `exclude`). VCF
text is not read by region: it is converted to BCF first.

### Personal sequences

`ConsensusReader` applies the records of a BCF, read by region through its
`.csi`, to a reference, as `bcftools consensus` 1.23 does, entirely in
memory. VCF text is converted to BCF first.

```python
from omicsweave.query import ConsensusOptions, ConsensusReader
from omicsweave.types import Haplotype, Samples

with ConsensusReader("cohort.bcf", genome, options=ConsensusOptions(include="QUAL>20")) as cohort:
    cohort.fetch(region, Haplotype("NA12878", "1")).text
    cohort.fetch(region.with_strand("-"), Haplotype("NA12878", "1"))  # reverse complement
    cohort.fetch(region, Samples(["NA12878", "NA12891"]))  # IUPAC codes
    batch = cohort.fetch_batch(region, [Haplotype(s, h) for s in cohort.samples for h in "12"])
```

Each region decodes only the records it takes, each once. A batch reads the
reference and the records of its region once, and returns every sequence
packed in one array. On a `-` region every sequence is the reverse
complement (IUPAC codes complemented); chains stay on the forward strand.

### Gene annotation

`omicsweave.annotation` reads GTF (GENCODE, Ensembl, RefSeq, UCSC) and GFF3,
plain, gzip or BGZF, whole and in parallel, into columnar tables.

```python
from omicsweave.annotation import Annotation, read_aliases

ann = Annotation.read("gencode.v44.annotation.gtf.gz")  # attributes="all" | [keys] | None
ann.genes.id, ann.genes.start, ann.transcripts.gene  # numpy / StringVector columns
ann.gene("TP53").transcripts                         # by id or name
ann.overlaps(region, level="exon")                   # rows of ann.exons
ann.report.problems                                  # orphans, strand conflicts, ...

coding = ann.subset(biotypes=["protein_coding"], tags=["basic"])
ann.introns(coding), ann.junctions(), ann.exon_bins(aggregate=True)
ann.promoters(2000, 200, chrom_lengths=lengths)      # IntervalSet, owner = transcript
ann.gene_bodies(upstream=2000)                       # Signac-style gene activity
ann.match_chroms(bam_reference_names)                # chr1 <-> 1, chrM <-> MT
ann.exonic_regions().overlaps(region, mode="within")  # also "overlaps", "containing"
ann.others.type.to_list()                            # GFF3 regions and other records
read_aliases("chromAlias.txt")                       # {"1": "chr1", ...} for aliases=
ann.write_gtf("copy.gtf")                            # byte-identical for GTF + all attributes
```

Genes and transcripts without records of their own (UCSC exports, GFF3
genes with only CDS children) are synthesized from their children.
Derived views are `IntervalSet`s with an `owner` column and their own
overlap index; `IntervalSet.from_arrays` builds one from peaks or bins.

## Writing

```python
from omicsweave.formats import BigBedSchema, BigBedWriter, BigWigWriter, FastaWriter, TwoBitWriter

with FastaWriter("out.fa.gz") as fasta:  # bgzip; saves .fai and .gzi
    fasta.add("chr1", "ACGT...")

with TwoBitWriter("out.2bit", ["chr1"]) as twobit:  # names first, then each in order
    twobit.add("chr1", "ACGTNNacgt")  # N and lowercase runs become blocks

with BigWigWriter("out.bw", {"chr1": 248_956_422}, threads=8) as bw:  # compress on 8 workers
    bw.add("chr1", starts, ends, values)  # columns or scalars

with BigBedWriter("out.bb", chroms, BigBedSchema(field_count=6)) as bb:
    bb.add("chr1", starts, ends, ["name\t0\t+", ...])
```

`HicWriter` writes binned maps and `HicBuilder` bins contact pairs at
every resolution, as Juicer's `pre` does; `blocks=HicBlocks(...)` sets the
bins per block. Writers compress on the calling thread unless given
`threads`; the file is the same either way. They finish on leaving the
`with` block and return the companion files they wrote;
`read_sidecar(path)` reads back what a sidecar records of a file.

## Converting

```sh
ow convert hg38.fa.gz hg38.2bit
ow convert cohort.vcf.gz cohort.bcf --index -t 16 -m 8G
ow convert reads.sam reads.bam --index
ow convert unsorted.bam sorted.bam --sort --index
ow convert reads.cram reads.bam --reference hg38.fa.gz --index
ow convert peaks.bed peaks.bb --chrom-sizes hg38.fa.fai
ow convert coverage.bedGraph coverage.bw --chrom-sizes hg38.chrom.sizes
ow convert phyloP.wig.gz phyloP.bw --chrom-sizes hg38.2bit
ow convert k562.mcool k562.hic --resolutions 1000,5000
ow convert k562.pairs.gz k562.hic -t 16 -m 16G
ow convert merged.txt k562.hic --chrom-sizes hg38.chrom.sizes --resolutions 5000,25000
ow convert pileup.bed calls.bb
ow convert pileup.bed 5mC.bw --mod m --value percent --strand +
ow convert genome.txt genome.out --from fasta --to bgzip-fasta
```

The source's format is read from its content and the target's from its
name; `--from` and `--to` state them. A BAM converts to a BAM only by
sorting it (`--sort`): the records, header and BGZF blocks are those
`samtools sort --no-PG` writes. A CRAM becomes the BAM htslib writes
from it (`samtools view -b --no-PG`), slices decoded in parallel;
`--reference` gives the FASTA or 2bit it was compressed against.
`--dry-run` shows the memory plan
without converting, and `ow convert --help` lists every conversion. The
target and its indexes appear only once the whole conversion succeeded.

bedGraph and wig (`fixedStep`, `variableStep` and bedGraph lines) become
bigWig as UCSC `bedGraphToBigWig` and `wigToBigWig` read them: values are
stored as given, input need not be sorted, and overlapping items are an
error. Without `--chrom-sizes`, each chromosome's length is its largest
item end.

4DN pairs and Juicer short text (`str1 chr1 pos1 frag1 str2 chr2 pos2
frag2`, with an optional score) become `.hic` binned as Juicer `pre` bins
them: every row counts once (or its score), positions are taken as written,
and nothing is filtered. The chromosomes come from the `#chromsize` lines or
`--chrom-sizes` (which short text needs); rows on other chromosomes are
skipped with a warning that counts them. The default resolutions are
`pre`'s (2.5 Mb down to 5 kb); coarser maps are summed from the finest, so
every resolution must be a multiple of it.

bedMethyl from modkit `pileup` (either delimiter) or ENCODE becomes a bigBed
that keeps every column, with an autoSql naming them, or a bigWig of one
column: `--value` is `percent` (default), `coverage`, `score`, or one of
modkit's counts (`nmod`, `ncanonical`, `nother_mod`, `ndelete`, `nfail`,
`ndiff`, `nnocall`), stored as written. `--mod` picks a modification code
(required when the file has several) and `--strand` a strand; selected items
that still overlap are an error naming the option to add.

The same from Python:

```python
report = ow.convert.run("cohort.vcf.gz", "cohort.bcf", index=True, threads=16, memory="8G")
report.target, report.companions  # Bcf('cohort.bcf'), (Csi('cohort.bcf.csi'),)
plan = ow.convert.plan("hg38.fa.gz", "hg38.2bit", memory="1G")  # every option checked
plan.memory, plan.budget  # nothing read yet
plan.run()
ow.convert.run("pileup.bed", "5hmC.bw", mod="h", value="coverage")
ow.convert.run("k562.pairs.gz", "k562.hic").warnings  # also issued as ConvertWarning
```

`threads` defaults to every core the process may run on (CPU affinity and
cgroup limits included), `memory` to 4 GiB.

Sorts that outgrow the memory budget spill to disk beside the target (or
to `tmp_dir`); the output is the same either way.

## Quantifying BAM

`ow quant` reads coordinate-sorted BAMs once and writes every product asked
for: coverage as bigWig, count matrices as 10x directories (with velocity
and labeling layers), single-cell ATAC fragments with their tabix index,
splice junctions, metabolic-labeling statistics, DNA methylation, Hi-C
contacts, allele counts at SNVs, and `summary.json`.

```sh
# single-cell RNA: genes x UMIs x CB (STARsolo Gene-like raw matrix), plus coverage
ow quant possorted.bam -o out --gtf genes.gtf --genes single-cell --coverage -t 16 -m 16G
# bulk RNA, featureCounts-style (genes x read pairs x BAM, counts.tsv), reverse-stranded
ow quant a.bam b.bam c.bam -o counts --gtf genes.gtf --genes bulk --pairs --library 2
# single-cell ATAC: fragments.tsv.gz + .tbi, metrics, peak / tile / gene-activity matrices
ow quant atac.bam -o atac --gtf genes.gtf --atac --peaks peaks.bb --tiles 5000
# stranded 5' end coverage, CPM, and ChIP bins x fragments
ow quant rna.bam -o tracks --coverage --coverage-unit five-prime --strands split --library 1 --scale cpm
ow quant chip.bam -o chip --chip-bins 1000 --matrix mtx
# splicing: junctions.tab (STAR SJ.out.tab columns, motifs from the genome), a
# junction x CB matrix, and spliced / unspliced / ambiguous layers (velocyto)
ow quant possorted.bam -o out --gtf genes.gtf --reference genome.fa.gz \
    --genes single-cell --velocity --junctions --junction-matrix single-cell
# SLAM-seq / 4sU: T>C per read (conversions.tsv.gz, mismatches.tsv) and a
# labeled layer (UMIs with >= 2 T>C), SNPs masked
ow quant slam.bam -o slam --gtf genes.gtf --reference genome.2bit --variants snps.bcf \
    --genes single-cell --labeled 2 --labeling --min-base-quality 20
# methylation (Bismark XM, MM/ML, or bwa-meth reads against the genome):
# cytosines.tsv.gz + .tbi, beta.bw, coverage.bw, mbias.tsv, 100 kb bins x CB
ow quant bs.bam -o meth --reference genome.fa.gz --methylation --contexts CpG,CHG \
    --merge-cpg --methylation-bins 100000 --methylation-rows CB
# Hi-C (after samtools fixmate -m): contacts.pairs.gz, contacts.hic, per-cell metrics
ow quant hic.bam -o hic --hic --hic-resolutions 1000000,100000,10000 --self-ligation 1000 \
    --hic-cells --genome hg38
# allele counts (cellsnp-lite style): SNV x CB ref/alt/other UMIs, variants.tsv.gz, and
# gene-level haplotype counts from a phased sample
ow quant possorted.bam -o snv --gtf genes.gtf --alleles snps.bcf --ase NA12878
```

Matrices go to `<product>/barcodes.tsv.gz`, `features.tsv.gz` and
`matrix.cpk` ([CellPack](https://github.com/CladeTeam/CellPack); `--matrix
mtx` or `csv` for `matrix.mtx.gz` / `matrix.csv.gz`). Peaks are read from
bigBed (`ow convert peaks.bed peaks.bb`). Further layers sit beside the
matrix (`spliced.cpk`, `unspliced.cpk`, `ambiguous.cpk`, `labeled.cpk`);
methylation, allele and haplotype matrices name every layer
(`methylated.cpk` / `total.cpk`, `ref.cpk` / `alt.cpk` / `other.cpk`,
`hap1.cpk` / `hap2.cpk` / `both.cpk`). The
reference is a bgzip FASTA with its `.fai` / `.gzi` (`ow convert genome.fa
genome.fa.gz`) or 2bit; without one, labeling and methylation read the `MD`
tags. Known
variants are a BCF with its CSI (`ow convert snps.vcf.gz snps.bcf --index`).
An unsorted BAM is refused with the command that sorts it; `--dry-run` shows
the memory plan. Results do not
depend on the thread count, and buffers that outgrow `-m` spill to disk.

The same from Python, with the annotation read once:

```python
from omicsweave import quant
from omicsweave.annotation import Annotation

task = quant.Task(
    ["possorted.bam"],
    "out",
    [
        quant.Coverage(strands="split", strand="forward"),
        quant.Counts.rna_single_cell(),  # presets: rna_bulk, atac, chip
        quant.Atac(peaks="peaks.bb"),
    ],
    annotation=Annotation.read("genes.gtf", attributes=None),
)
splicing = quant.Task(
    ["possorted.bam"],
    "splicing",
    [
        quant.Counts.rna_single_cell(velocity=True, labeled=quant.Labeled(min_conversions=2)),
        quant.Junctions(matrix=quant.JunctionMatrix.single_cell()),
        quant.Labeling(),
    ],
    annotation=Annotation.read("genes.gtf", attributes=None),
    reference="genome.fa.gz",
    variants="snps.bcf",
)
epigenome = quant.Task(
    ["bs.bam"],
    "meth",
    [quant.Methylation(matrices=[quant.MethylMatrix("bins", quant.Bins(100_000), rows="CB")])],
    reference="genome.fa.gz",
)
contacts = quant.Task(["hic.bam"], "hic", [quant.Contacts(barcode_tag="CB", min_distance=1000)])
snvs = quant.Task(
    ["possorted.bam"],
    "snv",
    [quant.Alleles("snps.bcf", genes=quant.AlleleGenes(sample="NA12878"))],
    annotation=Annotation.read("genes.gtf", attributes=None),
)
genes = Annotation.read("genes.gtf", attributes=None)
introns = genes.introns(genes.subset(biotypes=["protein_coding"]))
custom = quant.Counts.rna_bulk(quant.Regions(introns), name="introns")  # any IntervalSet
basic = quant.Counts.rna_bulk(subset=genes.subset(tags=["basic"]))  # genes of a subset
plan = quant.plan(task, threads=16, memory="16G")  # headers read, feature spaces built
plan.memory, plan.inputs  # fixed, pipeline, buffer; each BAM's references
report = plan.run()  # or quant.run(task, ...)
report.files, report.summary["products"]
```

Every product is a frozen dataclass whose defaults and presets come from
the Rust layer; the Rust layer also checks every option when the task is
planned.

Coverage of aligned bases matches `samtools depth`; bulk gene counts match
featureCounts (single-end and read pairs, every strandedness, `-O`); the
fragments and cytosine indexes are the ones tabix builds; the `.hic` file is
the one `HicBuilder` (and so `ow convert contacts.pairs.gz`) writes for the
same contacts. The definitions (Tn5 shift, cut sites, duplicate collapsing,
per-barcode TSS counts, junction overhangs and motifs, velocity classes, T>C
conversions, methylation calls and contexts, contact ends and filters,
allele units) and the differences from STARsolo, velocyto, featureCounts,
Bismark / MethylDackel, pairtools and cellsnp-lite are in
[the quantification guide](docs/quantification.md).

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
- Dense results run 5' to 3' of the key's strand; sparse results keep
  genomic order. `"."` and `"?"` read as `"+"`.
- Failures raise the exceptions of `omicsweave.errors`, each also derived
  from the closest built-in (`UnknownChromError` is a `LookupError`,
  `IoError` an `OSError`, ...).
- Readers release the GIL while decoding, so threads may share one; they
  pickle as the files and options they were opened with.
- Sidecars (`.sidecar`) describe a file's data blocks for loaders; reading
  by region does not need them. Build one with `build_sidecar()`, or
  `sidecar=True` on a writer or a conversion; `read_sidecar()` reads it.

## Fidelity

Every format is tested against the tool it replaces, built from source in
`refs/` and called in process: `samtools faidx` / `bgzip` indexes, htslib
BCF and BAM encoding and their `.csi` / `.bai` indexes, BCF, BAM and CRAM
region queries and `samtools index` (`.crai` included), CRAM to BAM, `bcftools consensus` 1.23, libBigWig, lib2bit and
straw. [The architecture](docs/design.md), [format semantics](docs/rust-api.md),
[conversion](docs/conversion.md) and [quantification](docs/quantification.md)
document the design and deliberate differences.

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
| `docs/` | [user guides, API semantics and development documentation](docs/README.md) |

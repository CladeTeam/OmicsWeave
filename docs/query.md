# 按区域查询

[格式转换](convert.md) · [定量分析](quant.md) · [安装](../README.md#安装)

给定一段基因组坐标，OmicsWeave 可以返回这段区域的序列、信号、区间、比对、变异或接触矩阵。查询直接返回包含 NumPy 数组的结果对象。

```python
from omicsweave.query import SequenceReader, TrackReader
from omicsweave.types import Region

region = Region("chr1", 1_000_000, 1_131_072)
with SequenceReader("genome.fa.gz") as genome, TrackReader("atac.bw") as track:
    sequence = genome.fetch(region)
    signal = track.signal(region)

print(sequence.text[:20])
print(sequence.bases.shape)  # (131072,)，uint8 ASCII 碱基
print(signal.values.shape)   # (131072,)，float32；没有数据的位置为 NaN
```

这里的 `genome.fa.gz` 必须是 bgzip FASTA，并已有 `.fai` 和 `.gzi`。文件准备方法见 [FASTA](formats/fasta.md)。示例中的路径、染色体和样本名需要替换为自己的数据。

## 基本流程

所有查询都遵循同样的四步：

```text
文件（路径或格式对象） ──► Reader ──► 查询 key ──► 结果对象（NumPy 数组）
 genome.fa.gz              SequenceReader   Region        Sequence.bases
 atac.bw                   TrackReader      Spliced       Signal.values / SignalBins.values
 contacts.hic              ContactReader    ContactWindow ContactMatrix.values
```

1. **准备文件。** 每类数据使用一种可随机访问的格式，外部索引预先建好。手头是文本格式时先用 [convert](convert.md) 转换一次。
2. **打开 reader。** reader 读取 header 和索引，之后的查询只读取需要的数据块。打开时会按内容核对格式，传错文件会立即报错。
3. **构造 key。** 连续区域用 `Region`（也可以直接写 tuple），拼接区间用 `Spliced`，二维接触窗口用 `ContactWindow`。所有坐标都是 0-based 半开。
4. **取出数组。** 结果对象保留原始 key，并以 NumPy 数组或 `StringVector` 提供数据，可以直接拼入 batch。

不确定该用哪个 reader 时，`query.open` 会按文件内容选择：

```python
from omicsweave import query

with query.open("signal.bw") as reader:  # 返回 TrackReader
    print(type(reader).__name__, reader.chroms)

query.open("peaks.bed")
# InvalidError: bed is not read by region; convert it to bigBed with
# `ow convert peaks.bed peaks.bb --chrom-sizes <genome.sizes>`
```

对不能直接查询的格式，错误信息会给出可以照抄的转换命令。

## 选择 reader

| 需要的数据 | Reader 与方法 | 支持的输入 | 查询 key | 返回结果 |
|---|---|---|---|---|
| 参考序列 | `SequenceReader.fetch` | 普通 FASTA、bgzip FASTA、2bit | `Region` / `Spliced` | `Sequence` |
| 逐碱基信号 | `TrackReader.signal` | bigWig | `Region` / `Spliced` | `Signal` |
| 分箱信号统计 | `TrackReader.bins` | bigWig | `Region` / `Spliced` + 分箱大小 | `SignalBins` |
| 原始信号区间 | `TrackReader.intervals` | bigWig | `Region` | `SignalIntervals` |
| 分链信号 | `StrandedTrackReader` | 正、负链各一个 bigWig | 带 `+` / `-` 的 key | 与 `TrackReader` 相同 |
| 重叠的区间记录 | `IntervalReader.fetch` | bigBed | `Region` | `Intervals` |
| 接触矩阵或非零像素 | `ContactReader.matrix` / `pixels` | `.hic` v8 / v9 | `ContactWindow` | `ContactMatrix` / `ContactPixels` |
| 比对记录 | `AlignmentReader.fetch` | BAM、CRAM 3.0 / 3.1 | `Region` + 可选字段 | `Alignments` |
| 基因型 | `VariantReader.genotypes` | BCF | `Region` + 可选样本 | `Genotypes` |
| 变异字段与过滤结果 | `VariantReader.fetch` | BCF | `Region` + 字段、过滤表达式 | `Variants` |
| 应用变异后的序列 | `ConsensusReader.fetch` / `fetch_batch` | BCF + 参考序列 | `Region` + 等位选择 | `ConsensusSequence` / `ConsensusBatch` |

也可以让 `omicsweave.query.open(path)` 根据文件内容选择普通 reader。分链信号和 consensus 需要额外输入，使用对应的 reader 创建。

以下格式需要先[转换](convert.md)，不能直接交给区域 reader：

| 手头的文件 | 准备用于查询的文件 | 格式说明 |
|---|---|---|
| 普通 gzip FASTA | bgzip FASTA + `.fai` + `.gzi`，或 2bit | [FASTA](formats/fasta.md)、[2bit](formats/2bit.md) |
| BED | bigBed | [BED](formats/bed.md)、[bigBed](formats/bigbed.md) |
| bedGraph、wig | bigWig | [bedGraph](formats/bedgraph.md)、[wig](formats/wig.md)、[bigWig](formats/bigwig.md) |
| bedMethyl | bigBed（所有字段）或 bigWig（一个数值字段） | [bedMethyl](formats/bedmethyl.md) |
| VCF，包括压缩 VCF | BCF + `.csi` | [VCF](formats/vcf.md)、[BCF](formats/bcf.md) |
| SAM，包括压缩 SAM | BAM + `.bai` / `.csi` | [SAM](formats/sam.md)、[BAM](formats/bam.md) |
| `.cool`、`.mcool` | `.hic` | [cool](formats/cool.md)、[mcool](formats/mcool.md) |
| 4DN pairs、Juicer short | `.hic` | [pairs](formats/pairs.md)、[Juicer short](formats/juicer-short.md) |

[GTF](formats/gtf.md) 和 [GFF3](formats/gff3.md) 通过 `Annotation.read` 载入，用于查找基因、转录本和生成查询 key，见下文。

## 查询 key

### 连续区域：`Region`

```python
from omicsweave.types import Region

region = Region("chr1", 100, 200)
reverse = region.with_strand("-")
print(len(region))  # 100

# 接受 Region 的查询也接受这些 tuple
key = ("chr1", 100, 200)
stranded_key = ("chr1", 100, 200, "-")
```

坐标从 0 开始，右端不包含在区间内：`[100, 200)` 覆盖 100 个位置。`start` 和 `end` 必须是非负整数，小于 `2**32`，且 `start <= end`。染色体名使用文件中实际保存的名称，例如 `chr1` 与 `1` 在普通 reader 中不会自动互换。

`strand` 默认 `"."`，也可以是 `"+"`、`"-"` 或 `"?"`。方向影响如下：

| 结果 | `strand="-"` 时的行为 |
|---|---|
| 序列、consensus | 反向互补，按 5′ → 3′ 返回 |
| 逐碱基信号、分箱统计 | 按负链方向排列 |
| 稠密接触矩阵 | 对应的行轴或列轴反转 |
| 区间、信号 runs、比对、变异、接触 pixels | 保留基因组坐标和原有记录顺序，不作链过滤 |

普通 reader 将 `"."` 和 `"?"` 按正向读取。`StrandedTrackReader` 要求明确的 `"+"` 或 `"-"`，因为它需要据此选择信号文件。

### 多段区域：`Spliced`

`Spliced` 把同一染色体上的多个区间拼成一段，适合外显子序列或转录本信号。

```python
from omicsweave.types import Spliced

transcript = Spliced("chr1", [(100, 150), (300, 350)], strand="-")
print(len(transcript))  # 100，不包含两个区间之间的内含子

with SequenceReader("genome.fa.gz") as genome:
    mrna = genome.fetch(transcript)
    print(mrna.text)  # 拼接后按负链方向返回

print(transcript.to_genome([0, 49, 50, 99]))  # [349, 300, 149, 100]
print(transcript.to_local([349, 100]))        # [0, 99]
```

区间必须按基因组坐标递增排列，彼此不重叠，可以首尾相接。长度是各段长度之和；负链也使用递增的输入区间。`to_local` 对不属于这些区间的位置返回 `-1`。

`SequenceReader.fetch`、信号的 `signal` / `bins` 支持 `Spliced`；区间、比对、变异和 consensus 查询使用连续 `Region`。

### 二维窗口：`ContactWindow`

```python
from omicsweave.types import ContactWindow

window = ContactWindow(
    ("chr1", 0, 1_000_000),
    ("chr2", 0, 500_000),
    resolution=5_000,
)
print(window.shape)  # (200, 100)
```

行、列可以来自不同染色体。窗口包含每个与区域重叠的 bin，因此区域边界不必对齐 bin。`resolution` 必须是文件中已有的分辨率，可用 `reader.resolutions()` 查询。默认 `unit="BP"`；文件含酶切片段分辨率时，也可使用 `unit="FRAG"`，此时区域坐标表示片段编号。

## 返回结果

结果对象保留原始 key。序列和信号使用 `.locus`，区域记录使用 `.region`，接触数据使用 `.window`。

| 结果类型 | 主要字段 | dtype 与形状 |
|---|---|---|
| `Sequence` | `bases`、`text` | `uint8 [长度]`；`text` 为 ASCII 字符串 |
| `Signal` | `values` | `float32 [长度]` |
| `SignalBins` | `values`、`stats`、`edges` | `float32 [统计量数, bin 数]`；edges 是局部坐标边界 |
| `SignalIntervals` | `start`、`end`、`value` | `uint32 [N]`、`uint32 [N]`、`float32 [N]` |
| `Intervals` | `start`、`end`、`annotation` | 两列 `uint32 [N]`；其余 BED 字段为 `StringVector` |
| `ContactMatrix` | `values` | `float32 [行 bin 数, 列 bin 数]` |
| `ContactPixels` | `row`、`col`、`value` | 两列 `uint32 [N]`；值为 `float32 [N]` |
| `Alignments` | 固定比对列及选取的可变列 | 见下文 |
| `Genotypes` | 位置、等位、`gt`、`phased` | `gt` 为 `int8 [记录数, 样本数, 倍性]` |
| `Variants` | `genotypes`、`qual`、`filters`、`info`、`format` | 按 BCF 声明的字段类型返回 |
| `ConsensusSequence` | `bases`、`text`、可选 `chain` | `uint8 [输出长度]`，indel 可改变长度 |
| `ConsensusBatch` | `bases`、`offsets`、`chains` | 多条变长序列打包；`batch[i]` 取一条 |

`StringVector` 用于变长字符串或字节列：`.to_list()` 转为 Python 列表，`.bytes(i)` 读取单条原始字节，`.to_arrow()` 转为 Arrow 字符串数组（需要安装 `pyarrow`）。比对质量和 aux 等二进制字段应按原始字节读取。

## 序列

```python
with SequenceReader("genome.2bit", soft_mask=True) as genome:
    print(genome.chroms)          # 染色体名 → 长度
    print(genome.length("chr1"))
    seq = genome.fetch(("chr1", 100, 200, "-"))
    print(seq.text)
```

FASTA 按原样保留大小写。2bit 默认返回大写；`soft_mask=True` 将软屏蔽区域返回为小写。序列查询必须在参考序列范围内，不会用 N 补齐越界部分。参考文件、索引和自定义索引路径见 [FASTA](formats/fasta.md) 与 [2bit](formats/2bit.md)。

## 信号与分箱统计

```python
with TrackReader("signal.bw") as track:
    dense = track.signal(("chr1", 0, 1000))
    bins = track.bins(("chr1", 0, 1000), 100, ("mean", "coverage"))
    print(bins.values.shape)  # (2, 10)
    print(bins["mean"])
    print(bins.edges)        # [0, 100, ..., 1000]
    runs = track.intervals(("chr1", 0, 1000))
    print(runs.start, runs.end, runs.value)
```

逐碱基信号中没有数据的位置为 NaN，超过已知染色体末端的部分同样为 NaN。`intervals` 返回文件中存储的 runs，并裁剪到查询区域。

`bins` 直接汇总区间数据，支持一次计算多个统计量：

| 统计量 | 含义 | 空 bin |
|---|---|---|
| `mean` | 有数据碱基的均值 | NaN |
| `min` / `max` | 有数据碱基的最小值 / 最大值 | NaN |
| `sum` | 按覆盖碱基数加权的总和 | NaN |
| `std` | 有数据碱基的样本标准差；只有一个碱基时为 0 | NaN |
| `coverage` | 有数据碱基占 bin 的比例 | 0 |

例如，100 bp 的 bin 中只有 20 bp 的值为 5，`mean=5`、`sum=100`、`coverage=0.2`。若把缺失碱基视为 0，该 bin 的均值为 `sum / 100 = 1`。

默认 `phase=0`，从 key 的局部位置 0 开始分箱，末尾可以是不完整 bin。`phase` 必须满足 `0 <= phase < resolution`，局部位置 `x` 属于 `(x + phase) // resolution`。正向连续区域若需对齐全基因组网格，可设置 `phase=region.start % resolution`。

默认 `exact=True` 汇总原始数据；`exact=False` 可利用文件中的 zoom 层加速大 bin 查询，结果是估计值，没有合适 zoom 时使用原始数据。`track.zoom_reductions` 列出已有 zoom；`track.summary` 是文件的整体汇总。

正负链分开存储时：

```python
from omicsweave.query import StrandedTrackReader

with StrandedTrackReader("rna.plus.bw", "rna.minus.bw") as rna:
    values = rna.signal(("chr1", 100, 200, "-")).values
```

负链 key 选择 minus 文件并反转位置顺序，数值的正负号按文件原样保留。

## 区间记录

```python
from omicsweave.query import IntervalReader

with IntervalReader("peaks.bb") as peaks:
    hits = peaks.fetch(("chr1", 100, 200))
    print(hits.start, hits.end)
    print(hits.annotation.to_list())
    print(peaks.schema)
```

返回与区域重叠的完整记录，按起点排列；记录的坐标不会裁剪。`annotation` 保存第 4 列起的字段，以 tab 分隔，不会自动拆成 name、score 等列。BED12 的 block 字段也保存在这里。格式和写入方法见 [bigBed](formats/bigbed.md)。

## 接触矩阵

```python
from omicsweave.query import ContactReader

with ContactReader("contacts.hic") as contacts:
    print(contacts.resolutions())
    window = ContactWindow(("chr1", 0, 1_000_000), ("chr1", 0, 1_000_000), 5_000)
    matrix = contacts.matrix(window)
    pixels = contacts.pixels(window)
    print(matrix.values.shape)  # (200, 200)
    print(pixels.row, pixels.col, pixels.value)
```

返回原始接触计数，不应用 KR、VC 或其他归一化。稠密矩阵没有接触的位置为 0；`pixels` 只返回非零位置，`row` / `col` 是 bin 起点的基因组坐标，不是窗口内的矩阵下标。文件信息见 [.hic](formats/hic.md)。

## 比对记录

```python
from omicsweave.query import AlignmentReader

with AlignmentReader("reads.bam") as reads:
    hits = reads.fetch(("chr1", 100, 200), fields=("name", "seq", "cigar"))
    print(hits.start, hits.end, hits.flag, hits.mapq)
    if len(hits):
        print(hits.names[0], hits.seqs[0], hits.cigars.string(0))
```

默认只解码固定列：

| 字段 | dtype | 含义 |
|---|---|---|
| `start` / `end` | `int64` | 比对起点与参考跨度的右端，0-based 半开 |
| `flag` | `uint16` | SAM FLAG |
| `mapq` | `uint8` | 比对质量 |
| `mate_tid` | `int32` | mate 的参考序列编号；`-1` 表示没有 |
| `mate_start` | `int32` | mate 的 0-based 起点；`-1` 表示没有 |
| `tlen` | `int32` | SAM TLEN |

`mate_tid` 对应 `reader.references` 中的顺序。可选字段名为 `name`、`seq`、`qual`、`cigar`、`aux`，或用 `fields="all"` 全部读取；未请求的字段为 `None`。质量是原始 Phred 字节，不是 FASTQ 的 ASCII 编码；缺失质量为 `0xff`。CIGAR 用 `.string(i)` 查看，aux 保存 BAM 二进制标签。

查询按文件顺序返回区域记录，不按 MAPQ、duplicate、secondary 等自动过滤，也不会因为负链 key 而反向互补读段。需要计数或覆盖度时使用 [quant](quant.md)。

CRAM 查询通常需要压缩时的参考：

```python
with AlignmentReader("reads.cram", reference="genome.fa.gz") as reads:
    hits = reads.fetch(("chr1", 100, 200), "all")
```

只读固定列、CIGAR、名称、质量时可以不提供参考；embedded-reference 或无参考压缩的 CRAM 也可自行恢复序列。校验和、MD/NM 和索引说明见 [CRAM](formats/cram.md)。

## 变异与基因型

```python
from omicsweave.query import VariantReader

with VariantReader("cohort.bcf") as variants:
    print(variants.samples)
    gt = variants.genotypes(("chr1", 0, 100_000), ["sample1", "sample2"])
    print(gt.start, gt.ref.to_list(), gt.alt.to_list())
    print(gt.gt.shape)  # (记录数, 2, 该区域的最大倍性)

    records = variants.fetch(
        ("chr1", 0, 100_000),
        samples=["sample1"],
        info=["DP", "AF"],
        format=["AD"],
        include="QUAL>30",
    )
    print(records.qual, records.info["AF"], records.format["AD"])
```

字段必须存在于输入 BCF 的 header。`samples` 默认选择所有样本，并按请求顺序返回。`genotypes` 返回位置、ID、REF、逗号连接的 ALT、`gt` 和 `phased`；`gt` 中 0 是 REF、1 是第一个 ALT，依此类推。`GT_MISSING=-1` 表示缺失等位，`GT_ABSENT=-2` 表示该位置没有这一倍性的等位或记录没有 GT。`phased` 是 `bool [记录数, 样本数]`。

`fetch` 额外返回 QUAL、FILTER 和指定的 INFO / FORMAT。`include` / `exclude` 使用 bcftools 表达式，二者互斥。数值 INFO 的形状是 `[记录数, width]`，FORMAT 是 `[记录数, 样本数, width]`；整数为 `int32`，浮点为 `float32`，INFO Flag 为 `bool [记录数]`。`width` 是此次查询中的最大字段长度。

变长字段用缺失值和向量结束值补齐：整数分别为 `INT_MISSING`、`INT_END`；浮点均显示为 NaN，可通过位模式 `FLOAT_MISSING` / `FLOAT_END` 区分。这些常量可从 `omicsweave.types` 导入。字符串返回 `StringVector`，FORMAT 字符串按 `记录下标 * 样本数 + 样本下标` 排列。文件准备见 [BCF](formats/bcf.md)。

## 个体化序列：consensus

```python
from omicsweave.query import ConsensusOptions, ConsensusReader
from omicsweave.types import AllAlt, Haplotype, Samples

with SequenceReader("genome.fa.gz") as genome:
    with ConsensusReader(
        "cohort.bcf", genome, options=ConsensusOptions(include="QUAL>20", chain=True)
    ) as cohort:
        region = Region("chr1", 1000, 2000)
        sequence = cohort.fetch(region, Haplotype("sample1", "1"))
        batch = cohort.fetch_batch(region, [Haplotype("sample1", "1"), Haplotype("sample1", "2")])
        print(sequence.text, sequence.chain)
        print(batch[0].text, batch[1].text)
```

| 等位选择 | 含义 |
|---|---|
| `AllAlt()`，也是默认值 | 对每条记录取第一个 ALT，不参考样本 GT |
| `Samples(["sample1", ...])` | 根据这些样本的 GT 生成 IUPAC 序列 |
| `Haplotype("sample1", "1")` / `"2"` | 取样本指定的单倍型等位 |

`Haplotype` 还支持 `R`、`A`、`LR`、`LA`、`SR`、`SA`、`I`、`<N>pIu` 等 bcftools 选择方式。`ConsensusOptions` 可以设置 `include` / `exclude`、缺失基因型字符 `missing`、无变异位置字符 `absent`、插入/缺失/SNV 标记、BED mask 和 `chain` 输出。

consensus 的区域是参考坐标；插入和缺失会改变结果长度。批量结果按 `offsets` 打包，不能假定是等长二维数组。默认 `regions_overlap="record"`，也可设为 `"start"` 或 `"variant"` 改变区域记录的重叠判定。

## 从注释生成 key

```python
from omicsweave.annotation import Annotation

annotation = Annotation.read("genes.gtf", attributes=None)
transcript = annotation.locus("ENST00000269305")  # exons，返回 Spliced
cds = annotation.locus("ENST00000269305", "cds")

with SequenceReader("genome.fa.gz") as genome:
    mrna = genome.fetch(transcript)
    coding_sequence = genome.fetch(cds)
```

`part` 支持 `exons`、`cds`、`utr5`、`utr3`；CDS 包含 stop codon。`annotation.loci()` 批量生成 key，`annotation.overlaps(region, level="gene")` 返回重叠基因的行号。注释也可生成 introns、promoters 等区间供 quant 使用。输入与导出方法见 [GTF](formats/gtf.md) 和 [GFF3](formats/gff3.md)。

## 结果存储、重复查询与并发

默认返回的数组由结果持有，后续查询、清空缓存或关闭 reader 不会使它失效。使用 `with` 可以在查询结束后关闭文件。

固定大小查询可选 `out=`，适合把结果直接写进已有批次数组：

```python
import numpy as np

buffer = np.empty(1000, dtype=np.float32)
with TrackReader("signal.bw") as track:
    result = track.signal(("chr1", 0, 1000), out=buffer)
    assert result.values is buffer
```

`out` 必须可写、C-contiguous，dtype 和形状完全匹配。序列使用 `uint8 [长度]`，signal 使用 `float32 [长度]`，bins 使用 `float32 [统计量数, bin 数]`，contact matrix 使用 `float32 [行数, 列数]`。修改或复用 buffer 会改变引用它的结果。

重复读取相邻区域时，可以给 reader 配置缓存：

```python
from omicsweave.cache import BlockCache

cache = BlockCache(max_bytes=256 * 1024**2)
with TrackReader("signal.bw", cache=cache) as track:
    first = track.signal(("chr1", 0, 1000))
    second = track.signal(("chr1", 500, 1500))
print(cache.stats())
```

bgzip FASTA、2bit、bigWig、bigBed、hic、BAM、CRAM、BCF 支持缓存；普通 FASTA 不接受非空缓存。缓存预算限制已接纳的数据，不是进程总内存上限。使用期间不要替换数据、索引或参考文件。

解码会释放 GIL。多个线程可共享 reader；需要独立解码状态时使用 `reader.clone()`，clone 共享文件和缓存。reader 可以 pickle，并在 DataLoader worker 中重新打开文件；缓存不会随之传递，`BlockCache` 必须在每个进程内新建，也不能跨 fork 使用。convert 和 quant 的 Plan 是进程内对象，不能 pickle。

## 查看文件信息

reader 打开后即可查看文件的元信息，不需要发起查询：

| Reader | 属性与方法 |
|---|---|
| 所有 reader | `source`（打开的格式对象）、`chroms`（名称 → 长度）、`closed` |
| `SequenceReader` | `length(chrom)` |
| `TrackReader` | `summary`（全文件统计）、`zoom_reductions` |
| `StrandedTrackReader` | `plus`、`minus`（两条链各自的 `TrackReader`） |
| `IntervalReader` | `schema`（字段数、标准字段数、autoSql） |
| `ContactReader` | `version`、`genome`、`attributes`、`maps`（有矩阵的染色体对）、`resolutions(unit)`、`sites(chrom)`、`fragment_count(chrom)` |
| `AlignmentReader` | `references`（按编号排列，对应 `mate_tid`）、`header_text`、`index` |
| `VariantReader` | `samples`、`index`；`chroms` 中 header 未给长度的 contig 为 `None` |
| `ConsensusReader` | `samples`、`reference`、`options`、`index` |

## 在 PyTorch DataLoader 中使用

reader 的设计目标就是放进 `Dataset`。推荐的写法是：Dataset 只保存格式对象（轻量、可 pickle），在每个 worker 第一次取样本时打开 reader 和缓存，之后所有查询复用它们。

```python
import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset

from omicsweave.cache import BlockCache
from omicsweave.formats import BgzipFasta, BigWig
from omicsweave.query import SequenceReader, TrackReader
from omicsweave.types import Region

LOOKUP = np.zeros((256, 4), dtype=np.float32)
for i, base in enumerate(b"ACGT"):
    LOOKUP[base, i] = LOOKUP[base + 32, i] = 1.0  # 大小写都编码


class Windows(Dataset):
    def __init__(self, genome: BgzipFasta, tracks: list[BigWig], regions: list[Region]):
        self.genome, self.tracks, self.regions = genome, tracks, regions
        self._readers = None  # 在 worker 中打开

    def _open(self):
        if self._readers is None:
            cache = BlockCache(max_bytes=256 << 20)  # 每个 worker 一个缓存
            self._readers = (
                SequenceReader(self.genome, cache=cache),
                [TrackReader(t, cache=cache) for t in self.tracks],
            )
        return self._readers

    def __len__(self):
        return len(self.regions)

    def __getitem__(self, i):
        genome, tracks = self._open()
        region = self.regions[i]
        x = LOOKUP[genome.fetch(region).bases]                     # [L, 4]
        y = np.stack([t.bins(region, 128).values[0] for t in tracks])  # [轨道数, L/128]
        return torch.from_numpy(x), torch.from_numpy(np.nan_to_num(y))


regions = [Region("chr1", s, s + 131_072) for s in range(1_000_000, 50_000_000, 65_536)]
loader = DataLoader(
    Windows(BgzipFasta("genome.fa.gz"), [BigWig("atac.bw"), BigWig("ctcf.bw")], regions),
    batch_size=32,
    num_workers=8,
    persistent_workers=True,          # reader 跨 epoch 保持打开
    multiprocessing_context="spawn",  # 原生 reader 和缓存在各 worker 内创建
)
```

要点：

- **格式对象可以 pickle，reader 也可以。** reader 被 pickle 时只保存打开参数，在 worker 中重新打开文件；`BlockCache` 不能 pickle，也不能跨 fork 使用，必须在每个进程内新建。
- **一个 worker 内共享缓存。** 序列和多条信号轨道共用一个 `BlockCache`，相邻或重叠的窗口复用已解码的数据块。
- **采样器在主进程。** 让 sampler 产出 `Region`，新 epoch 的窗口即可传给 persistent worker，无需重建 Dataset。
- **多线程场景**（例如在一个进程内用线程池取样）可以共享同一个 reader；需要并行解码时，每个线程使用 `reader.clone()`。
- **预分配 batch。** 固定大小的查询可以用 `out=` 直接写入 batch 数组的一行，省去一次拷贝。

完整示例见 [examples](../examples)：NonB 的轨道 / peak 加载器，以及按细胞组织多种测定、Hi-C tile 和表达量的 CENO 加载器。

## 错误类型

所有异常都继承自 `omicsweave.errors.OmicsWeaveError`，并同时继承对应的 Python 内置异常，可以按需捕获：

| 异常 | 也是 | 常见原因 |
|---|---|---|
| `UnknownChromError` | `LookupError` | 文件中没有该染色体名 |
| `InvalidRegionError` | `ValueError` | `start > end` |
| `RegionOutOfRangeError` | `IndexError` | 序列查询超出参考长度 |
| `OutputSizeError` | `ValueError` | `out=` 的形状或 dtype 不匹配 |
| `InvalidError` | `ValueError` | 参数或记录不合法；格式需要先转换 |
| `FormatError` | `ValueError` | 文本无法解析，附带 1-based 行号 |
| `IoError` | `OSError` | 文件或索引不存在、无法读取 |
| `BadMagicError`、`TruncatedError`、`InflateError` | — | 文件损坏、截断或不是预期格式 |
| `UnsupportedError` | `NotImplementedError` | 合法但未实现的输入，如 big-endian BBI、CRAM 2.x |
| `StaleSidecarError` | — | sidecar 与当前数据文件不匹配，需要重建 |

## 常见问题

| 现象 | 处理方法 |
|---|---|
| 提示缺少索引 | 显式调用对应格式的 `build_index()`；打开 reader 不会自动建索引 |
| 文件叫 `.fa.gz` 却不能查询 | 用 `formats.detect(path)` 确认是 bgzip FASTA；普通 gzip 需先转换 |
| 找不到染色体 | 检查 `reader.chroms`，使用文件中的实际名称 |
| 信号出现 NaN | 表示没有数据；按自己的分析需要决定是否填 0 |
| 接触窗口分辨率不可用 | 从 `reader.resolutions()` 选择已存储分辨率 |
| CRAM 序列无法恢复 | 提供压缩时对应的参考，并核对参考和索引 |
| 需要过滤比对或统计分子数 | 在返回列上自行筛选，或使用 [quant](quant.md) |

# 区域查询

[文档首页](README.md) · [格式与索引](formats.md) · [块缓存](query-cache.md)

## 选择 reader

```python
from omicsweave import Region
from omicsweave.query import SequenceReader, TrackReader

region = Region("chr1", 1000, 2024)
with SequenceReader("genome.fa.gz") as genome, TrackReader("signal.bw") as track:
    sequence = genome.fetch(region)         # sequence.bases: uint8 [1024]
    signal = track.signal(region)           # signal.values: float32 [1024]
    bins = track.bins(region, 128, ("mean", "coverage"))  # [2, 8]
    mean = bins["mean"]                    # float32 [8]
```

| Reader | 输入 | 方法 | 结果及主要列 |
|---|---|---|---|
| `SequenceReader` | FASTA / bgzip FASTA / 2bit | `fetch(locus)` | `Sequence.bases`：ASCII `uint8 [L]`；`.text` 转字符串 |
| `TrackReader` | bigWig | `signal(locus)` | `Signal.values`：`float32 [L]`，缺失为 NaN |
| `TrackReader` | bigWig | `bins(locus, resolution, stats)` | `SignalBins.values`：`float32 [S, B]` |
| `TrackReader` | bigWig | `intervals(locus)` | `start`、`end`、`value`；裁剪到查询跨度 |
| `StrandedTrackReader` | plus / minus 两个 bigWig | `signal`、`bins`、`intervals` | 按 key 的链选择轨道；无链 key 报错 |
| `IntervalReader` | bigBed | `fetch(region)` | `Intervals.start` / `end`：`uint32`，`annotation`：字符串列 |
| `ContactReader` | hic v8 / v9 | `matrix(window)` / `pixels(window)` | `ContactMatrix.values` 或 `row`、`col`、`value` |
| `AlignmentReader` | BAM / CRAM | `fetch(region, fields=...)` | `Alignments` 定长列与可选变长列 |
| `VariantReader` | BCF | `genotypes(region, samples=None)` / `fetch(region, samples=None, info=, format=, include=, exclude=)` | `Genotypes`；`Variants`（`Genotypes` 加 QUAL、FILTER、INFO / FORMAT 列） |
| `ConsensusReader` | BCF + 参考 | `fetch(region, alleles)` / `fetch_batch(region, alleles)` | `ConsensusSequence` / `ConsensusBatch` |

索引和版本要求见[格式与索引](formats.md)。`ow.open(path)` 根据内容选择单文件 reader。

## 坐标和方向

`Region(chrom, start, end, strand=".")` 使用 0-based、半开区间；坐标在 `[0, 2**32)` 中，`start <= end`。也可传 `(chrom, start, end)` 或含 strand 的四元 tuple。

`Spliced(chrom, spans, strand=".", name=None)` 描述同一染色体上升序、不重叠的多个跨度。序列和信号按跨度拼接；它不适用于 bigBed、alignment、variant、consensus 或 contact window 的区域参数。由注释构造转录本 key 见[基因注释](annotation.md)。

稠密序列、逐位信号、bins 按 key 的 5′→3′ 排列：负链翻转，序列还取互补。`"."` 和 `"?"` 按正链读取。稀疏结果保留基因组坐标，不因负链 key 自动过滤或翻转。

Sequence 跨出参考序列范围会报错；Track 超出染色体末端的部分缺失；bigBed 返回完整重叠特征而不裁剪。未知染色体和格式的边界差异见[Rust 格式语义](rust-api.md)。

## 信号汇总

支持 `mean`、`min`、`max`、`sum`、`std`、`coverage`，一次请求多个统计量共用解码。`mean` 的分母是有数据的碱基数；`coverage` 是覆盖比例。空 bin 除 coverage 为 0，其余为 NaN。

`phase=0` 时 bin 从拼接后的局部位置 0 开始。`phase` 移动 bin 网格，首尾可能为部分 bin。对正链连续区域，对齐基因组网格可用 `phase=region.start % resolution`。`exact=False` 允许用 zoom 近似，没有合适 zoom 时回退到全分辨率。数值兼容差异见[Rust 接口](rust-api.md)。

## Contact、alignment 与 genotype

```python
from omicsweave.query import AlignmentReader, ContactReader, VariantReader
from omicsweave.types import ContactWindow

with ContactReader("sample.hic") as contacts:
    window = ContactWindow(("chr1", 0, 100_000), ("chr1", 0, 100_000), 1000)
    matrix = contacts.matrix(window).values  # float32 [100, 100]，缺失为 0
    pixels = contacts.pixels(window)         # 坐标是 bin 起点，不是局部数组下标
with AlignmentReader("reads.bam") as reads:
    alignments = reads.fetch(region, fields=("name", "seq", "cigar"))
with VariantReader("cohort.bcf") as variants:
    genotypes = variants.genotypes(region, ["sample1"])
    records = variants.fetch(region, ["sample1"], info=["DP"], format=["AD"], include="QUAL>30")
```

Contact 仅返回原始计数，不读取 balance / expected；每条轴计入所有与区间重叠的 bin，负链轴翻转稠密矩阵。`ContactWindow.shape` 给出精确输出形状。FRAG 单位为酶切片段编号，需要文件中存在相应分辨率。

Alignment 默认只读固定列：`start`、`end`、`flag`、`mapq`、`mate_tid`、`mate_start`、`tlen`。`fields="all"` 或字段列表可读取 name、seq、qual、cigar、aux。查询不隐式过滤 duplicate、secondary 或低 MAPQ。CRAM 的固定字段、CIGAR、名字、质量无需参考；序列和生成 MD / NM 等路径可能需要 `reference=`。

GT 是 `int8 [records, samples, ploidy]`：非负值为等位编号，`-1` 缺失，`-2` 表示此倍性槽不存在；phase 为 `bool [records, samples]`。不额外计算 dosage 或 allele frequency。`fetch` 一次扫描同时解出 QUAL（缺失为 NaN）、FILTER 和所列 INFO / FORMAT 字段：数值按区间内最宽记录对齐为 `[records, width]`（FORMAT 为 `[records, samples, width]`），缺失值为 `INT_MISSING`（浮点 NaN），补位为 `INT_END`；`include` / `exclude` 是 bcftools 表达式。

## Consensus

```python
from omicsweave.query import ConsensusOptions, ConsensusReader, SequenceReader
from omicsweave.types import Haplotype

with SequenceReader("genome.fa.gz") as genome:
    with ConsensusReader("cohort.bcf", genome, options=ConsensusOptions(include="QUAL>20")) as cohort:
        first = cohort.fetch(region, Haplotype("sample1", "1"))
        batch = cohort.fetch_batch(region, [Haplotype("sample1", h) for h in ("1", "2")])
```

每次公开查询重新加载区域；同一次 batch 内共享参考与记录解码。Indel 可改变序列长度，输出不保证与 `len(region)` 相等。负链结果反向互补，chain 保持正链坐标。完整语义以移植的 bcftools consensus 规则为准，见[Rust 接口](rust-api.md)。

## 输出存储和生命周期

默认查询返回独立持有数据的结果；后续查询、缓存淘汰和关闭 reader 不使数组失效。字符串列使用 `StringVector`，可 `.to_list()` 或在安装 pyarrow 后 `.to_arrow()`；Arrow 视图共享底层存储。

固定大小的 `fetch`（sequence）、`signal`、`bins`、`matrix` 支持 `out=`：

```python
import numpy as np

with TrackReader("signal.bw") as track:
    storage = np.empty(len(region), dtype=np.float32)
    result = track.signal(region, out=storage)
    assert result.values is storage
```

`out` 必须可写、C 连续，shape 和 dtype 与查询完全匹配。复用 storage 会改变引用它的已有结果。缓存预算不包含这些数组和解码工作区。

## 线程、进程和异常

Reader 在解码期间释放 GIL；同一 reader 的可变工作区受锁保护，线程共享不代表同一上下文可以同时解码。需要独立执行状态时，缓存 reader 可 `clone()`。普通查询不要求创建缓存或先调用 describe / prefetch。

Reader 的 pickle 只保留重新打开的配置；缓存内容不随 pickle 传递。fork 后不要使用继承的 native reader，应在 DataLoader worker 内新建。最稳妥的用法是 dataset 只持有格式对象，在 worker 第一次请求时创建 reader 和缓存，见[示例](../examples/README.md)。

失败使用 `omicsweave.errors` 的异常，每类还继承对应的内置异常，例如 `IoError` 属于 `OSError`、`UnknownChromError` 属于 `LookupError`。已关闭 reader 的查询会报错，`close()` 可重复调用。

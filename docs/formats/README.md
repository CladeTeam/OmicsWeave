# 格式总览

[按区域查询](../query.md) · [格式转换](../convert.md) · [定量分析](../quant.md) · [项目首页](../../README.md)

基因组数据散落在几十种文件格式中：有的为人阅读而设计（FASTA、BED、VCF、SAM），有的为快速随机访问而设计（bigWig、BAM、BCF、.hic），还有一类是专门帮助定位数据的索引（FAI、BAI、CSI）。这一组文档逐个介绍 OmicsWeave 涉及的格式：它保存什么、内部如何组织、有哪些容易踩的坑，以及 OmicsWeave 怎样读、写和转换它。

如果只想快速找到“手里的文件该怎么用”，先看下面的分类表；需要理解格式本身时，再进入单篇文档。

## 按数据类型浏览

| 数据 | 文本 / 交换格式 | 可按区域查询的格式 | 所需索引 | OmicsWeave reader |
|---|---|---|---|---|
| 参考序列 | [FASTA](fasta.md)（普通 gzip） | [FASTA](fasta.md)（普通、bgzip）、[2bit](2bit.md) | [FAI](fai.md)、[GZI](gzi.md)；2bit 自带 | `SequenceReader` |
| 数值信号 | [bedGraph](bedgraph.md)、[wig](wig.md) | [bigWig](bigwig.md) | 内置 | `TrackReader`、`StrandedTrackReader` |
| 区间与 peaks | [BED](bed.md) | [bigBed](bigbed.md) | 内置 | `IntervalReader` |
| 甲基化调用 | [bedMethyl](bedmethyl.md) | bigBed（全部字段）、bigWig（一个数值列） | 内置 | `IntervalReader`、`TrackReader` |
| 读段比对 | [SAM](sam.md) | [BAM](bam.md)、[CRAM](cram.md) | [BAI](bai.md) / [CSI](csi.md)、[CRAI](crai.md) | `AlignmentReader` |
| 变异与基因型 | [VCF](vcf.md) | [BCF](bcf.md) | [CSI](csi.md) | `VariantReader`、`ConsensusReader` |
| 染色质接触 | [cool](cool.md)、[mcool](mcool.md)、[4DN pairs](pairs.md)、[Juicer short](juicer-short.md) | [.hic](hic.md) | 内置 | `ContactReader` |
| 基因注释 | [GTF](gtf.md)、[GFF3](gff3.md) | 载入内存后查询 | 不需要 | `Annotation` |
| 加载元信息 | — | [Sidecar](sidecar.md)（bigWig / bigBed / .hic 的伴随文件） | — | `read_sidecar` |

OmicsWeave 的原则是：**每类数据选一个能高效随机访问的格式作为查询格式，其他格式通过一次转换进入这个格式。** 这样 reader 只需要做好少数几种二进制格式，转换在数据准备阶段完成一次，训练时不再反复解析文本。

```text
bedGraph ─┐
wig ──────┼──► bigWig ──► TrackReader
bedMethyl ┘

BED ──────┬──► bigBed ──► IntervalReader
bedMethyl ┘

VCF ─────────► BCF + CSI ──► VariantReader / ConsensusReader
SAM ─────────► BAM + BAI/CSI ──► AlignmentReader ◄── CRAM + CRAI
CRAM ────────► BAM（需要时）
cool / mcool / pairs / Juicer short ──► .hic ──► ContactReader
普通 gzip FASTA ──► bgzip FASTA + FAI + GZI，或 2bit ──► SequenceReader

坐标排序 BAM ──► ow quant ──► bigWig / 计数矩阵 / .hic / pairs / ...
```

转换的具体组合和选项见 [格式转换](../convert.md)，从 BAM 计算覆盖度、计数矩阵等产物见 [定量分析](../quant.md)。

## 支持矩阵

| 格式 | 识别 | 区域查询 | 转换输入 | 转换输出 | Python 直接写入 | 建索引 |
|---|---|---|---|---|---|---|
| [FASTA](fasta.md) 普通 / bgzip | 是 | 是 | 是 | bgzip | `FastaWriter` | `build_index()` |
| FASTA 普通 gzip | 是 | — | 是 | — | — | — |
| [2bit](2bit.md) | 是 | 是 | — | 是 | `TwoBitWriter` | 自带 |
| [bigWig](bigwig.md) | 是 | 是 | — | 是 | `BigWigWriter` | 自带 |
| [bigBed](bigbed.md) | 是 | 是 | — | 是 | `BigBedWriter` | 自带 |
| [.hic](hic.md) v8 / v9 | 是 | 是 | — | 是 | `HicWriter`、`HicBuilder` | 自带 |
| [BAM](bam.md) | 是 | 是 | 是（排序） | 是 | — | `build_index()` |
| [CRAM](cram.md) 3.0 / 3.1 | 是 | 是 | 是 | — | — | `build_index()` |
| [BCF](bcf.md) | 是 | 是 | — | 是 | — | `build_index()` |
| [VCF](vcf.md)、[SAM](sam.md) | 是 | — | 是 | — | — | — |
| [BED](bed.md)、[bedGraph](bedgraph.md)、[wig](wig.md)、[bedMethyl](bedmethyl.md) | 是 | — | 是 | — | — | — |
| [cool](cool.md)、[mcool](mcool.md) | 是 | — | 是 | — | — | — |
| [4DN pairs](pairs.md)、[Juicer short](juicer-short.md) | 是 | — | 是 | quant 输出 pairs | — | — |
| [GTF](gtf.md)、[GFF3](gff3.md) | 由 `Annotation.read` 判断（`detect` 会当作 BED） | 载入后查询 | — | — | `write_gtf`、`write_bed12` | 不需要 |
| [FAI](fai.md)、[GZI](gzi.md) | 无 magic，不识别 | — | — | 随 FASTA 生成 | 随 `FastaWriter` 生成 | — |
| [BAI](bai.md)、[CSI](csi.md)、[CRAI](crai.md) | 是 | — | — | 随转换生成 | — | — |
| [Sidecar](sidecar.md) | 是 | — | — | 可选生成 | writer 的 `sidecar=True` | `build_sidecar()` |

“识别”指 `omicsweave.formats.detect(path)` 能从文件内容判断格式；它不看扩展名。

## 几个贯穿所有格式的概念

### 坐标系

不同格式记录位置的方式不同，这是最常见的错误来源：

| 约定 | 含义 | 使用它的格式 |
|---|---|---|
| 0-based 半开 `[start, end)` | 第一个碱基是 0，右端不包含 | BED、bedGraph、bigWig、bigBed、bedMethyl、BAM 二进制、2bit |
| 1-based 闭区间 `[start, end]` | 第一个碱基是 1，两端都包含 | GTF、GFF3、VCF 的 POS、SAM 的 POS、wig 的 start |
| 1-based 单点位置 | 一个碱基的位置 | 4DN pairs、Juicer short 的 pos |

**OmicsWeave 的所有查询 key 和返回结果都使用 0-based 半开坐标。** 文本中的 1-based 位置在读取、转换时统一换算：GTF 中 `start=101, end=200` 载入后是 `[100, 200)`；VCF `POS=101` 的 SNV 返回 `start=100`。

```text
碱基:        A   C   G   T   A
0-based:     0   1   2   3   4      BED "chr1 1 3" = CG
1-based:     1   2   3   4   5      GTF "2 3"      = CG
```

### 压缩：普通 gzip 与 BGZF

`.gz` 后缀可能代表两种完全不同的压缩：

- **普通 gzip**：整个文件是一个连续的压缩流，想读中间的内容必须从头解压。
- **BGZF**（blocked gzip，由 `bgzip` 生成）：文件由许多独立的 gzip 块组成，每块解压后不超过 64 KiB，并在块头记录自身大小。任意一块都能单独解压，配合索引即可随机访问。

BGZF 仍是合法的 gzip，`zcat` 可以直接读，所以从扩展名和普通工具上看不出区别。BAM、BCF 内部都使用 BGZF；FASTA 只有 bgzip 版本能按区域查询。用 `formats.detect` 可以区分：

```python
from omicsweave.formats import detect

print(detect("genome.fa.gz"))  # BgzipFasta(...) 或 GzipFasta(...)
```

文本类转换输入（FASTA、BED、bedGraph、wig、bedMethyl、VCF、SAM、pairs、Juicer short）都接受纯文本、普通 gzip 和 BGZF 三种存储。

### 索引与伴随文件

随机访问需要知道“目标数据在文件的哪个位置”。有两种做法：

- **内置索引**：bigWig、bigBed、.hic、2bit 把索引写在文件内部，单个文件即可查询。
- **外部索引**：FASTA、BAM、BCF、CRAM 把索引放在旁边的独立文件中（`.fai`、`.gzi`、`.bai`、`.csi`、`.crai`），习惯上命名为“数据文件名 + 后缀”。

外部索引与数据文件的字节位置绑定：数据文件被改写、重新压缩或排序后，旧索引立即失效，必须重建。OmicsWeave 的 reader 默认在数据文件旁查找索引，也可以用参数显式指定；**打开 reader 永远不会自动建索引**，需要时显式调用格式对象的 `build_index()`，或在转换时加 `--index`。

### 排序

区域索引依赖记录按坐标排序。BAM、BCF、CRAM 建索引前必须已经按坐标排序；BED、bedGraph、wig 转换时 OmicsWeave 会自动排序；VCF → BCF、SAM → BAM 保留输入顺序，不排序。未排序的 BAM 可以用 `ow convert in.bam out.bam --sort` 排序。

### 染色体名称

`chr1` 和 `1`、`chrM` 和 `MT` 在不同来源的数据中混用很常见。普通 reader 只接受文件中实际保存的名称，不做自动换算，可以先查看 `reader.chroms`。注释与 BAM 在 quant 中匹配时会自动处理 `chr1` / `1`、`chrM` / `MT`，并支持额外的别名表，见 [GTF](gtf.md)。

## 格式对象

`omicsweave.formats` 中每种格式（以及每种存储方式）都有一个类，例如 `Fasta`、`BgzipFasta`、`GzipFasta`。格式对象只是“这个路径上是某种格式的文件”的声明：

```python
from omicsweave.formats import BigWig, Bam, detect, FORMATS

bw = BigWig("signal.bw")       # 不打开文件，也不检查内容
print(bw.path, bw.exists(), bw.format)

bam = detect("reads.bam")      # 按内容识别，返回 Bam(...)
index = bam.build_index()      # 需要文件操作的方法才真正读写文件

print(sorted(FORMATS))         # 全部格式名，例如 'bgzip-fasta'、'bigwig'
```

格式对象可以代替路径传给 reader、`convert.run` 和 writer；当自动识别不符合预期时（例如四列 BED 被识别为 bedGraph），用格式对象或 CLI 的 `--from` / `--to` 显式指定。格式对象可以 pickle，适合放进 DataLoader 的 Dataset 中，在每个 worker 里再打开 reader。

## 单篇文档

- 序列：[FASTA](fasta.md) · [2bit](2bit.md) · [FAI](fai.md) · [GZI](gzi.md)
- 信号与区间：[bigWig](bigwig.md) · [bedGraph](bedgraph.md) · [wig](wig.md) · [bigBed](bigbed.md) · [BED](bed.md) · [bedMethyl](bedmethyl.md)
- 比对：[SAM](sam.md) · [BAM](bam.md) · [CRAM](cram.md) · [BAI](bai.md) · [CRAI](crai.md)
- 变异：[VCF](vcf.md) · [BCF](bcf.md) · [CSI](csi.md)
- 接触：[.hic](hic.md) · [cool](cool.md) · [mcool](mcool.md) · [4DN pairs](pairs.md) · [Juicer short](juicer-short.md)
- 注释：[GTF](gtf.md) · [GFF3](gff3.md)
- OmicsWeave 伴随文件：[Sidecar](sidecar.md)

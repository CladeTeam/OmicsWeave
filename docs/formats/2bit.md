# 2bit

[格式总览](README.md) · [查询](../query.md) · [转换](../convert.md) · [FASTA](fasta.md)

2bit 是 UCSC 设计的**紧凑二进制参考序列格式**：每个碱基只占 2 bit，四个碱基装进一个字节；N 区域和小写软屏蔽区域不逐个存储，而是记成若干“起点 + 长度”的块。文件自带序列目录，不需要任何外部索引就能按区域读取。

在 OmicsWeave 中，2bit 与 [FASTA](fasta.md) 一样可交给 `SequenceReader` 查询，也可以作为 CRAM 解码、consensus、quant 的参考，以及转换时的 `chrom_sizes`。它体积约为未压缩 FASTA 的四分之一，查询只需一次读取，适合训练时大量随机读取参考序列。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 若干条命名的 DNA 序列 |
| 常见后缀 | `.2bit` |
| 存储 | 二进制，每碱基 2 bit；N 和软屏蔽记为区间块 |
| 字符集 | 仅 A、C、G、T、N，外加大小写（软屏蔽）信息 |
| 坐标 | 0-based 半开 |
| 随机访问 | 支持，文件内置目录 |
| 外部索引 | 不需要 |
| 典型来源 | UCSC 发布的 `hg38.2bit` 等；`faToTwoBit` |
| OmicsWeave | 区域查询、直接写入、FASTA → 2bit 转换、作为参考和染色体长度来源 |

## 格式由来

2bit 由 Jim Kent 为 UCSC 基因组浏览器设计，配套工具 `faToTwoBit`、`twoBitToFa`、`twoBitInfo` 随 UCSC Kent 工具发布。它的目标是让浏览器和 BLAT 等工具快速取出任意一段参考序列：DNA 只有四种碱基，2 bit 正好够用，而 N 与软屏蔽在真实基因组中都是成片出现的，用区间描述远比逐碱基标记省空间。

## 文件结构

```text
┌──────────────────────────────────────┐
│ header（16 字节）                     │  magic 0x1A412743、版本、序列数、保留字段
├──────────────────────────────────────┤
│ 目录：每条序列一项                     │  名称长度（1 字节）、名称、记录偏移
├──────────────────────────────────────┤
│ 序列记录 1                            │
│   dnaSize          序列长度           │
│   nBlockCount      N 块数             │
│   nBlockStarts[]   N 块起点           │
│   nBlockSizes[]    N 块长度           │
│   maskBlockCount   软屏蔽块数          │
│   maskBlockStarts[]                   │
│   maskBlockSizes[]                    │
│   reserved                            │
│   packedDna        ⌈dnaSize / 4⌉ 字节 │
├──────────────────────────────────────┤
│ 序列记录 2 ...                        │
└──────────────────────────────────────┘
```

- **header**：版本 0 的目录偏移为 32 bit，每条序列记录的起点必须不超过 `2^32 - 1`；最后一条记录可以延伸到 4 GiB 之外。版本 1 使用 64 bit 偏移。
- **目录**：序列名最长 255 字节，偏移指向该序列记录的起点。因为目录位于所有记录之前，写文件时必须预先知道全部序列名。
- **packedDna**：每字节 4 个碱基，高位在前，编码为 `T=00`、`C=01`、`A=10`、`G=11`。N 位置在 packedDna 中也占位（按 `T` 填充），真实含义由 N 块覆盖。

读取一个区域的过程：在目录中找到序列 → 读取记录头（N 块、屏蔽块）→ 计算区域在 packedDna 中的字节范围并一次读出 → 解码为碱基 → 需要时用屏蔽块转为小写 → 最后用 N 块覆盖为大写 `N`。

例如，序列 `ACGTNacgtN`：

| 内容 | 存储方式 |
|---|---|
| `ACGT` `acgt` 的碱基 | packedDna 中的 2-bit 码 |
| 第 4、9 位的 `N` | N 块 `(4, 1)`、`(9, 1)` |
| 第 5–8 位小写 | 软屏蔽块 `(5, 4)` |

## 关键概念

### 只保留 ACGTN

2bit 无法保存 IUPAC 模糊碱基。转换或写入时，A、C、G、T（不分大小写）以外的字符一律存成 N，与 UCSC 工具的行为一致。人类参考基因组中模糊碱基很少，但如果下游依赖这些字符（例如某些物种的参考或人工构建的序列），应保留 [FASTA](fasta.md)。

### 软屏蔽

软屏蔽是“哪些区域是重复序列”的标记，存储在屏蔽块中，与碱基本身分离。OmicsWeave 默认返回全大写碱基（相当于 `twoBitToFa -noMask`）；`soft_mask=True` 时屏蔽区域返回小写（相当于 `twoBitToFa`）。N 块优先于屏蔽块：即使某段 N 落在屏蔽区内，也始终返回大写 `N`。

### 序列长度即染色体长度

2bit 的记录头直接给出每条序列的长度，所以它也是一份可靠的染色体长度表。OmicsWeave 的转换选项 `chrom_sizes` / `--chrom-sizes` 可以直接接受 `.2bit`。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 格式识别 | `detect` 按 magic 识别为 `TwoBit` |
| 区域查询 | `SequenceReader.fetch`，返回 `uint8` ASCII 的 `Sequence` |
| 多段 key | 支持 `Spliced` |
| 链方向 | 负链返回反向互补 |
| 软屏蔽 | 默认大写；`soft_mask=True` 返回小写屏蔽区 |
| 版本 | 读取 version 0 和 version 1（64-bit 偏移） |
| 缓存与并发 | 可接 `BlockCache` 缓存解码后的固定长度片段；解码释放 GIL；可 `clone()`、pickle |
| 写入 | `TwoBitWriter`，可选 `wide_offsets=True` 写 version 1 |
| 转换 | 纯文本、gzip、bgzip FASTA → 2bit；文件需要 64-bit 偏移时自动写 version 1 |
| 作为参考 | CRAM 解码、`ConsensusReader`、quant 的 `reference`；转换的 `chrom_sizes` |

限制：

- 只读 little-endian 文件；big-endian 2bit 报 `BigEndianError`（常见工具生成的都是 little-endian）。
- 区域必须位于序列长度之内，越界报错，不会补 N。
- 当前 convert 不提供 2bit → FASTA。

## 查询

```python
from omicsweave.query import SequenceReader

with SequenceReader("genome.2bit") as genome:
    print(genome.chroms)                       # 名称 → 长度，按文件顺序
    seq = genome.fetch(("chr1", 100, 200))
    print(seq.text)                            # 全大写

with SequenceReader("genome.2bit", soft_mask=True) as genome:
    seq = genome.fetch(("chr1", 100, 200, "-"))
    print(seq.text)                            # 屏蔽区为小写，负链反向互补
```

查询接口与 FASTA 完全相同：`Region`、`Spliced`、负链、`out=`、`BlockCache` 都可以使用，见 [FASTA](fasta.md) 和[查询文档](../query.md)。同一段代码可以在 FASTA 和 2bit 之间切换，只是 2bit 不会返回 IUPAC 字符。

## 从 FASTA 转换

```sh
ow convert genome.fa genome.2bit
ow convert genome.fa.gz genome.2bit -t 4 -m 2G
```

```python
from omicsweave import convert

convert.run("genome.fa.gz", "genome.2bit", threads=4, memory="2G")
```

普通、gzip、bgzip FASTA 均可作为输入：ACGT 以外的字符变为 N，小写保留为软屏蔽块。因为目录必须写在所有记录之前，转换先获得全部序列名：输入旁边有可用的 `.fai` 时直接读取，否则多扫描一遍输入。版本选择以序列记录起点是否能用 32 bit 表示为准，并计入 N / 屏蔽块的空间；仅凭现有 `.fai` 的长度无法判断时，会额外流式扫描输入。每条序列的碱基在 N / 屏蔽块确定之前暂存，超出内存预算时写入临时文件；指定 `tmp_dir` 时使用该目录，否则使用输出文件所在目录（裸文件名使用当前目录）。

## 从内存写入

```python
from omicsweave.formats import TwoBitWriter

with TwoBitWriter("example.2bit", names=["chr1", "chr2"]) as writer:
    print(writer.next_name)          # 'chr1'
    writer.add("chr1", "ACGTNNacgt")
    writer.add("chr2", "TTGGCCAA")
```

- `names` 预先声明全部序列，`add` 必须按这个顺序逐条调用；`next_name` 给出下一条应添加的名称，全部添加后为 `None`。
- 序列可以是 `str`、`bytes` 或 `uint8` 数组。ACGT 以外的字符存为 N，N 和小写的连续区段自动记为块。
- 任何序列记录的起点可能超过 `2^32 - 1` 时，使用 `wide_offsets=True` 写出 version 1；writer 不会自动升级版本。版本 0 只检查每条记录的起点，因此最后一条记录可以跨过 4 GiB 边界。
- 离开 `with` 时补写目录并完成文件；块内发生异常时放弃未完成的文件。

## 2bit 与 FASTA 的选择

| | 2bit | bgzip FASTA |
|---|---|---|
| 体积 | 约为未压缩 FASTA 的 1/4 | 通常更小（gzip 压缩），取决于压缩等级 |
| 查询开销 | 一次读取，位运算解码 | 读取并解压若干 BGZF 块 |
| 字符集 | ACGTN + 软屏蔽 | 原样保存任意字符 |
| 需要的文件 | 单个 `.2bit` | `.fa.gz` + `.fai` + `.gzi` |
| 其他工具兼容 | UCSC 工具、py2bit 等；samtools 不读取 | samtools、bcftools 等通用 |

两者在 OmicsWeave 中的查询结果，除 IUPAC 字符和默认大小写外完全相同。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 结果全是大写 | 2bit 默认不返回软屏蔽；使用 `soft_mask=True` |
| 模糊碱基变成了 N | 2bit 无法保存 IUPAC 码；需要时使用 FASTA |
| 打开时报 `BigEndianError` | 文件是 big-endian 2bit，用 UCSC 工具转回 FASTA 再转换 |
| 写入时报名称顺序错误 | `add` 必须按 `names` 顺序调用，查看 `writer.next_name` |
| 序列记录起点超过 32-bit 范围，写入失败 | `TwoBitWriter` 设置 `wide_offsets=True`；FASTA 转换会自动选择版本 |

## 相关格式

- [FASTA](fasta.md)：2bit 的文本来源，保留完整字符集。
- [CRAM](cram.md)：可用 2bit 作为解码参考。
- [bigWig](bigwig.md)、[bigBed](bigbed.md)、[.hic](hic.md)：转换时可用 2bit 提供染色体长度。

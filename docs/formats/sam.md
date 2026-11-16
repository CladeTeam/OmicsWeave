# SAM

[格式总览](README.md) · [转换](../convert.md) · [BAM](bam.md) · [CRAM](cram.md) · [查询](../query.md)

SAM（Sequence Alignment/Map）是保存**读段比对结果**的文本格式：每行一条记录，说明一条测序读段比对到参考基因组的哪个位置、以什么方式比对、质量如何。比对软件（BWA、Bowtie2、STAR、minimap2 等）通常直接输出 SAM，它也是 [BAM](bam.md) 和 [CRAM](cram.md) 共同的数据模型：三者保存的是同一套字段，只是编码方式不同。

在 OmicsWeave 中，SAM 是**转换源**：先转成 BAM，再用 `AlignmentReader` 按区域读取，或交给 [quant](../quant.md) 计算覆盖度和计数矩阵。理解 SAM 的字段，也就理解了 BAM / CRAM 查询结果中每一列的含义。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 读段比对记录：位置、FLAG、MAPQ、CIGAR、mate 信息、序列、质量、可选标签 |
| 常见后缀 | `.sam`、`.sam.gz` |
| 存储 | tab 分隔文本；可用 gzip / bgzip 压缩 |
| 坐标 | POS / PNEXT 从 1 开始；OmicsWeave 结果统一为 0-based 半开 |
| 随机访问 | 不支持，需要从头读 |
| 外部索引 | 无 |
| 典型来源 | 比对软件的直接输出、`samtools view` |
| OmicsWeave | 格式识别、SAM → BAM 转换（可同时建索引）；不直接查询，不写出 SAM |

## 格式由来

SAM 规范由 1000 Genomes 项目在 2009 年提出（Li et al., *Bioinformatics* 2009），随 samtools / htslib 一起发展，现由 GA4GH 的 [hts-specs](https://github.com/samtools/hts-specs) 维护。它统一了此前各比对软件互不兼容的输出，使下游工具只需理解一种格式。

文本 SAM 便于查看和调试，但体积大、无法按区域读取，所以实际存储和分析通常使用二进制的 BAM 或参考压缩的 CRAM。

## 文件结构

一个 SAM 文件由 header 和比对记录两部分组成：

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:248956422
@SQ	SN:chr2	LN:242193529
@RG	ID:lib1	SM:sample1	PL:ILLUMINA
@PG	ID:bwa	PN:bwa	VN:0.7.17
read1	99	chr1	10001	60	10M	=	10151	160	ACGTACGTAC	IIIIIIIIII	NM:i:0	RG:Z:lib1
read1	147	chr1	10151	60	10M	=	10001	-160	TTGCATGCAA	IIIIIIIIII	NM:i:1	RG:Z:lib1
```

### Header

以 `@` 开头的行是 header，每行由两字母类型和若干 `TAG:VALUE` 组成：

| 行类型 | 作用 | 常见字段 |
|---|---|---|
| `@HD` | 文件级信息，必须在第一行 | `VN` 格式版本，`SO` 排序方式（`coordinate`、`queryname`、`unsorted`） |
| `@SQ` | 参考序列字典，顺序决定参考序列编号 | `SN` 名称，`LN` 长度，`AN` 别名，`M5` 序列 MD5 |
| `@RG` | read group：测序文库 / 样本 / lane | `ID`，`SM` 样本，`LB` 文库，`PL` 平台 |
| `@PG` | 处理过该文件的程序 | `ID`，`PN`，`VN`，`CL` 命令行 |
| `@CO` | 自由文本注释 | — |

`@SQ` 特别重要：记录中的 RNAME / RNEXT 都必须是这里声明的名称，BAM 中参考序列以它们在 `@SQ` 中的顺序编号（tid）。

### 11 个必需列

| 列 | 名称 | 含义 | 缺失值 |
|---|---|---|---|
| 1 | QNAME | 读段名；同一模板（配对的两端）同名 | `*` |
| 2 | FLAG | 按位组合的状态标志，见下表 | — |
| 3 | RNAME | 比对到的参考序列名 | `*` |
| 4 | POS | 最左比对碱基的位置，**1-based** | `0` |
| 5 | MAPQ | 比对质量，`-10·log10(比对位置错误的概率)` | `255` 表示不可用 |
| 6 | CIGAR | 读段与参考如何对齐，见下表 | `*` |
| 7 | RNEXT | mate 的参考序列名；`=` 表示与 RNAME 相同 | `*` |
| 8 | PNEXT | mate 的位置，1-based | `0` |
| 9 | TLEN | 模板（插入片段）长度，最左端为正、最右端为负 | `0` |
| 10 | SEQ | 读段序列（反向比对时为参考正链方向） | `*` |
| 11 | QUAL | 碱基质量，Phred+33 编码的 ASCII | `*` |

第 12 列起是可选标签，格式为 `TAG:TYPE:VALUE`。

### FLAG

| 位 | 十进制 | `SamFlag` | 含义 |
|---|---|---|---|
| `0x1` | 1 | `PAIRED` | 测序为配对模板 |
| `0x2` | 2 | `PROPER_PAIR` | 两端按比对软件的标准正确配对 |
| `0x4` | 4 | `UNMAPPED` | 本端未比对 |
| `0x8` | 8 | `MATE_UNMAPPED` | mate 未比对 |
| `0x10` | 16 | `REVERSE` | 本端比对到负链（SEQ 已反向互补） |
| `0x20` | 32 | `MATE_REVERSE` | mate 比对到负链 |
| `0x40` | 64 | `READ1` | 模板中的第一个读段 |
| `0x80` | 128 | `READ2` | 模板中的第二个读段 |
| `0x100` | 256 | `SECONDARY` | 次要比对（多重比对中非主要的一条） |
| `0x200` | 512 | `QC_FAIL` | 未通过质控 |
| `0x400` | 1024 | `DUPLICATE` | PCR 或光学重复 |
| `0x800` | 2048 | `SUPPLEMENTARY` | 补充比对（嵌合比对的其他片段） |

例如 `99 = 1 + 2 + 32 + 64`：配对、正确配对、mate 在负链、read 1。`omicsweave.types.SamFlag` 提供这些位的枚举，可用于在查询结果上筛选：

```python
from omicsweave.types import SamFlag

keep = (hits.flag & (SamFlag.UNMAPPED | SamFlag.SECONDARY | SamFlag.DUPLICATE)) == 0
```

### CIGAR

CIGAR 由“长度 + 操作”重复组成，例如 `3S5M2D4M`：

| 操作 | 含义 | 消耗读段 | 消耗参考 |
|---|---|---|---|
| `M` | 比对上（可能匹配或错配） | 是 | 是 |
| `I` | 相对参考的插入 | 是 | 否 |
| `D` | 相对参考的缺失 | 否 | 是 |
| `N` | 跳过参考区域，RNA-seq 中表示内含子 | 否 | 是 |
| `S` | 软剪切：碱基保留在 SEQ 中但未比对 | 是 | 否 |
| `H` | 硬剪切：碱基不在 SEQ 中 | 否 | 否 |
| `P` | padding，用于多序列比对 | 否 | 否 |
| `=` | 序列匹配 | 是 | 是 |
| `X` | 序列错配 | 是 | 是 |

比对在参考上的跨度是所有“消耗参考”操作的长度之和。`3S5M2D4M` 从 POS 开始覆盖 5 + 2 + 4 = 11 个参考碱基。

### 常见可选标签

| 标签 | 类型 | 含义 |
|---|---|---|
| `NM` | `i` | 与参考的编辑距离 |
| `MD` | `Z` | 错配和缺失位置的字符串，配合 SEQ 可以还原参考 |
| `AS` / `XS` | `i` | 比对分数 / 次优比对分数 |
| `NH` | `i` | 该读段的比对位置数，RNA-seq 中用于判断多重比对 |
| `RG` | `Z` | 所属 read group |
| `MC` / `MQ` | `Z` / `i` | mate 的 CIGAR / MAPQ（samtools fixmate 补齐） |
| `SA` | `Z` | 嵌合比对的其他片段 |
| `CB` / `UB` | `Z` | 单细胞条形码 / UMI（10x 等，校正后） |
| `XM` / `XG` | `Z` | Bismark 甲基化调用 / 转换方向 |
| `MM` / `ML` | `Z` / `B` | 碱基修饰及其概率（如 5mC） |

类型字母：`A` 单字符、`i` 整数、`f` 浮点、`Z` 字符串、`H` 十六进制、`B` 数值数组。

## 关键概念

### 1-based 与 0-based

SAM 文本中的 POS / PNEXT 从 1 开始，区间两端都包含。OmicsWeave 的所有查询结果都使用 0-based 半开坐标：SAM `POS=10001`、CIGAR `10M` 的记录，在 `AlignmentReader` 结果中是 `start=10000`、`end=10010`。

### 排序

SAM 本身可以是任意顺序。比对软件通常按输入读段顺序（同名读段相邻）输出；按区域查询和建索引则要求**按坐标排序**：先按 `@SQ` 中的参考顺序，再按 POS。`@HD SO:` 只是声明，不保证文件真的有序。

### 未比对读段

FLAG 含 `0x4` 的读段仍可能带 RNAME / POS：配对中未比对的一端通常放在其 mate 的位置，排序后与 mate 相邻。RNAME 为 `*` 的读段没有位置，排在坐标排序文件的末尾，按区域查询取不到它们。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 格式识别 | `formats.detect` 根据 header 行（`@HD`、`@SQ` 等）识别为 `Sam`；压缩文件为 `GzipSam` |
| 压缩输入 | 纯文本、普通 gzip、BGZF 都可以 |
| SAM → BAM | `ow convert reads.sam reads.bam`，并行编码，输出与 `samtools view -b` 写出的 BAM 数据流一致 |
| 同时建索引 | `--index` / `index=True` 写 BAI 或 CSI，要求输入已按坐标排序 |
| 区域查询 | 不直接支持，先转 BAM |
| quant | 不直接接受，先转 BAM |
| 写出 SAM | 不提供 BAM → SAM 或 SAM writer |

与 htslib 的刻意差异：

- RNAME 或 RNEXT 不在 `@SQ`（或其 `AN` 别名）中时报错。htslib 会警告并把记录改成未定位，静默丢掉名称。
- 含 NUL 字节的行报错。
- 每行末尾的一个 `\r` 总会去掉，便于处理 Windows 换行。

## 转为 BAM

```sh
# 已按坐标排序的 SAM：转换并同时建索引
ow convert sorted.sam sorted.bam --index

# 未排序 SAM：先转二进制，再排序并建索引
ow convert reads.sam reads.bam
ow convert reads.bam sorted.bam --sort --index
```

```python
from omicsweave import convert

convert.run("reads.sam.gz", "reads.bam", threads=4)
report = convert.run("reads.bam", "sorted.bam", sort=True, index=True, memory="2G")
print(report.target, report.companions)  # Bam('sorted.bam') (Bai('sorted.bam.bai'),)
```

转换保留输入记录顺序，不自动排序；`--sort` 只适用于 BAM → BAM。`--index` 默认在每条参考序列都不超过 `2**29` bp 时写 [BAI](bai.md)，否则写 [CSI](csi.md)；也可用 `--index bai` / `--index csi` 指定。`--level` 设置 BGZF 压缩等级（0–12）。

转换后的 BAM 可以用 `AlignmentReader` 按区域读取，详见 [BAM](bam.md)。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 报参考序列不在 header 中 | RNAME / RNEXT 未在 `@SQ` 声明；补全 header 后再转换 |
| `--index` 报未排序 | 先不带 `--index` 转 BAM，再 `--sort --index` |
| 没有 header 的 SAM 识别失败 | 识别依赖 header 行，且需要 `@SQ` 声明参考；先补上 header |
| 想得到过滤后的读段 | 查询结果不自动过滤；在返回的 `flag`、`mapq` 列上筛选，或使用 quant 的过滤规则 |
| 需要把 BAM 导出为 SAM 查看 | OmicsWeave 不提供导出，可用 `samtools view` |

## 相关格式

- [BAM](bam.md)：SAM 的 BGZF 二进制形式，OmicsWeave 查询和 quant 的输入。
- [CRAM](cram.md)：基于参考压缩的比对格式，可转为 BAM。
- [BAI](bai.md)、[CSI](csi.md)：坐标排序 BAM 的区域索引。

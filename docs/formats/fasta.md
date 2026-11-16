# FASTA

[格式总览](README.md) · [查询](../query.md) · [转换](../convert.md) · [2bit](2bit.md) · [FAI](fai.md) · [GZI](gzi.md)

FASTA 是保存**命名序列**的文本格式：参考基因组的每条染色体、转录组的每条转录本、蛋白质组的每条蛋白，都可以写成“一行名称 + 若干行序列”。它几乎是所有序列分析工具的通用输入，也是参考基因组最常见的发布形式。

在 OmicsWeave 中，FASTA 是参考序列的主要来源：`SequenceReader` 按区域读取碱基，返回 `uint8` ASCII 数组；CRAM 解码、consensus、quant 的剪接 motif / 甲基化 / 代谢标记也使用它作为参考。普通 FASTA 和 bgzip FASTA 借助 [FAI](fai.md)（以及 [GZI](gzi.md)）可直接随机访问；普通 gzip FASTA 需要先转换。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 若干条命名序列（DNA、RNA 或蛋白） |
| 常见后缀 | `.fa`、`.fasta`、`.fna`、`.fa.gz` |
| 存储 | 纯文本；或 bgzip（BGZF）压缩；或普通 gzip 压缩 |
| 坐标 | 文件本身没有坐标；查询使用 0-based 半开区间 |
| 随机访问 | 纯文本与 bgzip 借助索引支持；普通 gzip 不支持 |
| 外部索引 | 纯文本需要 `.fai`；bgzip 需要 `.fai` + `.gzi` |
| 典型来源 | UCSC、Ensembl、GENCODE、NCBI 发布的参考基因组 |
| OmicsWeave | 区域查询、建索引、直接写入、转换为 bgzip FASTA 或 2bit、作为其他功能的参考 |

## 格式由来

FASTA 来自 1985 年 Lipman 与 Pearson 的 FASTP / FASTA 序列比对程序，它的输入格式因为简单而被广泛沿用，成为事实标准。格式从未有过严格的正式规范，不同工具对细节（名称规则、换行、允许的字符）略有差异；如今通行的“可索引 FASTA”约定来自 samtools `faidx`：每条序列的各行长度一致，这样可以通过算术直接定位任意碱基。

## 文件结构

```text
>chr1 AC:CM000663.2 Homo sapiens chromosome 1
ACGTNacgt
RYK
>chr2
TTGGCCAA
```

| 部分 | 说明 |
|---|---|
| header 行 | 以 `>` 开头。`>` 之后到第一个空白字符是**序列名**（这里是 `chr1`），其后的描述被忽略 |
| 序列行 | header 之后、下一个 `>` 之前的所有行，拼接起来就是这条序列 |
| 换行 | 一条序列可以写成多行；可索引的 FASTA 要求除最后一行外各行长度相同 |

序列字符的常见含义：

| 字符 | 含义 |
|---|---|
| `A` `C` `G` `T` | 四种碱基 |
| `N` | 未知碱基，常见于端粒、着丝粒和 gap |
| `R` `Y` `K` `M` `S` `W` `B` `D` `H` `V` | IUPAC 模糊碱基，例如 `R` = A 或 G |
| 小写 `acgt` | **软屏蔽**（soft-mask）：通常表示 RepeatMasker 等工具标出的重复序列，碱基本身仍然有效 |

### 三种存储方式

| 文件对象 | 存储 | 区域查询 | 所需索引 |
|---|---|---|---|
| `Fasta` | 纯文本，`.fa` / `.fasta` / `.fna` | 支持 | `.fai` |
| `BgzipFasta` | BGZF（bgzip），常见 `.fa.gz` | 支持 | `.fai` 和 `.gzi` |
| `GzipFasta` | 普通 gzip，常见 `.fa.gz` | 先转换 | 无法随机访问 |

BGZF 把文本切成不超过 64 KiB 的独立 gzip 块，依次拼接。整体仍然是合法的 gzip 文件，`zcat` 照常可读；但因为每块可以单独解压，配合索引就能只解压需要的块。普通 gzip 是一个连续的压缩流，读取第 N 个字节必须从头解压，所以不能随机访问。

两者的扩展名都是 `.gz`，**不能凭扩展名区分**。用 `omicsweave.formats.detect(path)` 按内容判断：

```python
from omicsweave.formats import detect

print(detect("genome.fa.gz"))  # BgzipFasta('genome.fa.gz') 或 GzipFasta('genome.fa.gz')
```

## 关键概念

### 索引如何定位碱基

samtools 风格的随机访问依赖两个事实：每条序列在文件中的起始字节已知，且每行的碱基数和字节数固定。于是第 `pos` 个碱基（0-based）位于

```text
offset + (pos // 每行碱基数) * 每行字节数 + pos % 每行碱基数
```

这些数字就保存在 [FAI](fai.md) 中。bgzip 文件里 FAI 给出的是**解压后**文本的位置，还需要 [GZI](gzi.md) 把它换算到压缩块。因此：

- 行长不一致的 FASTA 无法建索引，需要先解析序列，再用 `FastaWriter` 重新排版；`ow convert` 转成 bgzip FASTA 时保留原始文本，不修复行宽。
- 修改序列或换行方式后，旧索引立即失效，必须重建。

### 序列名

序列名只取 header 中 `>` 之后的第一个词。区域查询使用的就是这个名字：文件写的是 `chr1`，查询就必须用 `chr1`，不会与 `1` 自动互换。同名序列在建索引时只保留第一条。

### 大小写与软屏蔽

FASTA 中的小写字母通常是软屏蔽标记，不代表不同的碱基。OmicsWeave 读取 FASTA 时按原样返回大小写；需要统一大写时由调用方处理：

```python
bases = genome.fetch(region).bases
upper = bases & 0xDF  # ASCII 字母转大写（仅适用于字母）
```

### 链方向

序列只按正链存储。负链 key 返回反向互补，按 5′ → 3′ 排列：IUPAC 码映射到互补集合的码（如 `R` ↔ `Y`），`N`、`S`、`W` 互补后不变，大小写保留，`U` 按 `T` 互补。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 格式识别 | `detect` 按内容区分 `Fasta`、`BgzipFasta`、`GzipFasta`（以 `>` 开头的文本，压缩方式由首字节判断） |
| 区域查询 | `SequenceReader.fetch`，返回 `Sequence`：`bases` 为 `uint8` ASCII，`text` 为字符串 |
| 多段 key | 支持 `Spliced`，可直接拼接外显子得到成熟转录本序列 |
| 链方向 | 负链返回反向互补 |
| 建索引 | `Fasta.build_index()` 写 `.fai`；`BgzipFasta.build_index()` 一次扫描写 `.fai` 和 `.gzi`，与 samtools / bgzip 生成的文件字节一致 |
| 写入 | `FastaWriter`：纯文本或 bgzip，边写边生成索引 |
| 转换 | 任意 FASTA（纯文本、gzip、BGZF）→ bgzip FASTA 或 2bit |
| 缓存与并发 | bgzip FASTA 可接 `BlockCache` 缓存解压后的块；解码释放 GIL；可 `clone()`、pickle |
| 作为参考 | CRAM 解码、`ConsensusReader`、quant 的 `reference` |

限制：

- reader 不会自动建索引，缺少索引时报 `IoError`。
- 区域必须位于序列长度之内，越界报 `RegionOutOfRangeError`，不会用 N 补齐；需要补齐时自行截断区域再拼接。
- 普通 gzip FASTA 不能直接查询，也不能通过补一个 `.gzi` 变成可查询文件。
- 纯文本 FASTA 不接受非空的 `BlockCache`（它直接读取文件，没有需要缓存的解压块）。

## 准备文件

### 已有可索引的文件：只建索引

```python
from omicsweave.formats import BgzipFasta, Fasta

Fasta("genome.fa").build_index()             # genome.fa.fai
BgzipFasta("genome.fa.gz").build_index()     # genome.fa.gz.fai 和 genome.fa.gz.gzi
```

按自己的文件类型选择其中一条。建索引会检查内容确实是对应的格式；行长不一致等问题会报错。

### 普通 gzip 或想重新压缩：转换

```sh
ow convert genome.fa genome.fa.gz              # 纯文本 → bgzip，自动写 .fai 和 .gzi
ow convert ordinary-gzip.fa.gz genome.bgz.fa.gz
ow convert genome.fa.gz genome.2bit            # 任意 FASTA → 2bit
```

```python
from omicsweave import convert

report = convert.run("ordinary-gzip.fa.gz", "genome.bgz.fa.gz", threads=4)
print(report.target, report.companions)  # BgzipFasta(...)，(Fai(...), Gzi(...))
```

转换到 bgzip FASTA 时原样压缩解压后的 FASTA 文本，不重排 header、碱基或换行；行宽不齐仍会在生成 FAI 时报告 `different line length`。输出与 `bgzip` 在同一压缩等级下的结果逐字节一致，且不随线程数改变；`--level` / `level=` 设置 libdeflate 压缩等级（0–12）。选择 bgzip FASTA 还是 2bit，见下文“选择存储方式”。

## 查询

```python
from omicsweave.query import SequenceReader

with SequenceReader("genome.fa.gz") as genome:
    print(genome.chroms)               # mappingproxy({'chr1': 248956422, ...})，按文件顺序
    print(genome.length("chr1"))
    seq = genome.fetch(("chr1", 0, 100))
    print(seq.text)                    # 'NNNN...'
    print(seq.bases.dtype, seq.bases.shape)  # uint8 (100,)
```

索引默认是数据文件旁边的 `<fasta>.fai` 和 `<fasta>.gzi`。放在其他位置时显式传入：

```python
SequenceReader("genome.fa.gz", fai="index/genome.fai", gzi="index/genome.gzi")
```

### 负链与拼接区间

```python
from omicsweave.types import Spliced

with SequenceReader("genome.fa.gz") as genome:
    minus = genome.fetch(("chr1", 1000, 1100, "-"))   # 反向互补
    mrna = genome.fetch(
        Spliced("chr1", [(1000, 1200), (1500, 1700)], strand="-")
    )                                                 # 400 nt，拼接后按负链 5′ → 3′
```

`Spliced` 也可由注释生成：`Annotation.locus("ENST...")` 返回转录本外显子的 key，见 [GTF](gtf.md)。

### 转成模型输入

`bases` 是 ASCII 码，可以用查找表直接变成整数编码或 one-hot：

```python
import numpy as np

lut = np.full(256, 4, dtype=np.int8)          # 其他字符 → 4
for i, b in enumerate(b"ACGT"):
    lut[b] = i
    lut[b | 0x20] = i                          # 小写同样映射
tokens = lut[seq.bases]                        # int8 [长度]
one_hot = np.eye(5, dtype=np.float32)[tokens][:, :4]
```

### 批量与 DataLoader

固定长度窗口可用 `out=` 直接写进批次数组，避免额外拷贝：

```python
batch = np.empty((len(regions), 1000), dtype=np.uint8)
with SequenceReader("genome.fa.gz") as genome:
    for i, region in enumerate(regions):
        genome.fetch(region, out=batch[i])
```

reader 可 pickle，在 DataLoader worker 中会重新打开文件；滑动窗口等重复读取相邻区域的场景可以给 bgzip FASTA 配置 `BlockCache`。详见[查询文档](../query.md)。

## 从内存写入

```python
from omicsweave.formats import FastaWriter

with FastaWriter("example.fa.gz", line_width=60) as writer:
    writer.add("chr1", "ACGTNacgt")
    writer.add("chr2", b"TTGGCCAA")
# 离开 with 后得到 example.fa.gz、example.fa.gz.fai、example.fa.gz.gzi
```

- `compress` 默认按路径是否以 `.gz` 结尾决定，也可以显式指定 `True` / `False`；压缩时写 bgzip，`level` 为 libdeflate 等级 0–12，默认 6。
- `line_width` 是每行碱基数，默认 60；`0` 表示每条序列写成一行。
- 序列可以是 `str`、`bytes` 或 `uint8` 数组；名称必须唯一，序列非空，碱基限于 `0x21`–`0x7E` 的可打印非空白 ASCII，且不能包含 `>`。
- 索引根据写出的布局直接生成，不再回读文件；`close()` 返回写出的伴随文件 `(Fai, Gzi)` 或 `(Fai,)`。
- 块内发生异常时放弃未完成的文件。

## 选择存储方式

| 需求 | 推荐 |
|---|---|
| 需要保留 IUPAC 模糊碱基或原始字符 | bgzip FASTA |
| 需要与 samtools、bcftools 等工具共用参考 | bgzip FASTA（或纯文本 FASTA） |
| 只关心 ACGTN 和软屏蔽，想要更小、更快的文件 | [2bit](2bit.md) |
| 大文件长期存放 | bgzip FASTA 或 2bit，不要用普通 gzip |

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| `.fa.gz` 打不开，提示是 `GzipFasta` | 这是普通 gzip，先 `ow convert` 成 bgzip FASTA 或 2bit |
| 提示找不到 `.fai` | reader 不自动建索引；调用 `build_index()` 或指定 `fai=` |
| 建索引报行长不一致 | 同一序列除末行外行宽不齐；解析序列后用 `FastaWriter` 统一行宽，`ow convert` 不修复此问题 |
| 找不到染色体 | 名称取 header 第一个词；查看 `genome.chroms` |
| 区域越界报错 | 查询不会补 N；先按 `genome.length(chrom)` 截断 |
| 结果中有小写 | 文件中的软屏蔽；按需要转大写 |
| 换了 FASTA 后结果异常 | 旧的 `.fai` / `.gzi` 与新文件不匹配，重建索引 |

## 相关格式

- [FAI](fai.md)、[GZI](gzi.md)：FASTA 随机访问所需的索引。
- [2bit](2bit.md)：更紧凑的参考序列格式，可由 FASTA 转换得到。
- [CRAM](cram.md)：解码读段序列时需要与压缩时一致的参考 FASTA。
- [BCF](bcf.md)：`ConsensusReader` 把变异应用到 FASTA 参考上生成个体序列。

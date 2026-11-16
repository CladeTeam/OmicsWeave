# GZI：BGZF 块偏移索引

[格式总览](README.md) · [FASTA](fasta.md) · [FAI](fai.md) · [查询](../query.md)

`.gzi` 是 bgzip 定义的二进制索引，记录 BGZF 文件中每个压缩块的起点在**压缩文件**和**解压后文本**中各自的位置。有了它，就能把“解压后第 N 个字节”换算成“从哪个压缩块开始解压”。

在 OmicsWeave 中，GZI 与 [FAI](fai.md) 一起让 bgzip FASTA 支持区域查询：FAI 把碱基坐标换算成解压后文本的字节范围，GZI 再把这段范围换算成需要读取并解压的压缩块。

## 速览

| 项目 | 说明 |
|---|---|
| 作用 | 解压后偏移 → BGZF 压缩块位置 |
| 命名 | `<file>.gzi`，例如 `genome.fa.gz.gzi` |
| 存储 | 二进制，little-endian 64 位整数 |
| 配套文件 | bgzip FASTA 的 `.fai` |
| 生成工具 | `bgzip -i` / `bgzip -r`、`samtools faidx`、`BgzipFasta.build_index()`、`FastaWriter`、FASTA → bgzip 转换 |
| OmicsWeave | 建索引、查询 bgzip FASTA 时读取 |

## 背景：BGZF

BGZF（Blocked GNU Zip Format）是 htslib 使用的分块 gzip：数据被切成解压后不超过 64 KiB 的片段，每段单独压缩成一个完整的 gzip member，依次拼接，文件末尾是一个空的 EOF 块。每个块的 gzip header 中有一个额外字段（`BC`）记录该块的压缩长度。

整个文件仍是合法的 gzip，`zcat`、`gunzip` 可以照常解压；但因为每块可以单独解压，只要知道块的位置，就能跳到文件中间读取。BAM、BCF、tabix 索引的 VCF 也都使用 BGZF。普通 gzip 是单个连续的压缩流，没有块边界，无论补什么索引都不能随机访问。

## 文件结构

```text
┌────────────────────────────────────────┐
│ count        u64   条目数               │
├────────────────────────────────────────┤
│ compressed_0 u64   │ uncompressed_0 u64 │  第 2 个块的起点
│ compressed_1 u64   │ uncompressed_1 u64 │  第 3 个块的起点
│ ...                                    │
└────────────────────────────────────────┘
```

每个条目是一对位置：该块在压缩文件中的字节偏移，以及它的第一个字节在解压后文本中的偏移。起始 `(0, 0)` 隐含、不写入文件。OmicsWeave 建索引时只为非空块记录位置，首个非空块对应隐含的 `(0, 0)`；空块（包括 EOF）不写条目。读取索引时要求压缩偏移严格递增、解压偏移不减，允许空块对应相同的解压偏移；相邻条目的解压偏移之差不超过 64 KiB。

## 原理：GZI 如何配合 FAI

读取 bgzip FASTA 的 `chr1:[start, end)` 分三步：

1. **FAI**：用序列的 OFFSET 和行布局算出这段碱基在解压后文本中的字节范围 `[a, b)`。
2. **GZI**：二分查找，找到解压后起点 ≤ `a` 的最后一个块，以及覆盖到 `b` 的最后一个块，得到需要读取的压缩字节范围。
3. **解压**：读出这些块并逐块解压，从第一个块中跳过 `a − 块起点` 个字节，取出文本并去掉换行。

```text
碱基坐标 ──FAI──▶ 解压后字节 [a, b) ──GZI──▶ 压缩块 k … m ──inflate──▶ 碱基
```

因为块是独立的，相邻查询命中的块可以复用。OmicsWeave 给每个块一个稳定编号；显式传入 `BlockCache` 时可缓存解压后的块，减少滑动窗口等重复查询的解压工作。默认 `cache=None`，不启用共享块缓存。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 建索引 | `BgzipFasta.build_index()` 一次扫描，完整解压并校验各块以生成 `.fai`，同时生成 `.gzi`；Rust `Gzi::build` 单独建 GZI 时只读块头和 `ISIZE` 尾字段，不解压 |
| 随写随建 | 压缩的 `FastaWriter`、FASTA → bgzip 转换根据写出的块直接生成 |
| 查询 | `SequenceReader(path, gzi=...)`；默认使用 `<fasta>.gzi` |
| 路径对象 | `Gzi.beside(path)` 返回 `<path>.gzi` 的格式对象，不执行建索引 |
| 兼容性 | 读取 bgzip / samtools 生成的 GZI，包括省略空块条目的情况；生成的 GZI 与 bgzip 一致 |

限制：

- `.gzi` 没有固定签名，也没有专门的探测分支；文本探测兜底可能将其识别为 `Bed`，应使用 `Gzi(path)` 显式声明索引类型。
- OmicsWeave 只在 bgzip FASTA 查询中使用 GZI。BAM、BCF 通过 [BAI](bai.md) / [CSI](csi.md) 中的虚拟偏移定位，不需要 GZI。

## 使用

### 创建

```python
from omicsweave.formats import BgzipFasta

fai, gzi = BgzipFasta("genome.fa.gz").build_index()
print(fai.path, gzi.path)  # genome.fa.gz.fai genome.fa.gz.gzi
```

`ow convert genome.fa genome.fa.gz` 和 `FastaWriter("genome.fa.gz")` 会自动生成 `.fai` 和 `.gzi`。

### 查询时指定

```python
from omicsweave.query import SequenceReader

with SequenceReader(
    "genome.fa.gz", fai="index/genome.fai", gzi="index/genome.gzi"
) as genome:
    sequence = genome.fetch(("chr1", 0, 100))
```

默认路径为 `<fasta>.gzi`。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| bgzip FASTA 提示缺少 `.gzi` | 调用 `BgzipFasta(...).build_index()`，或用 `gzi=` 指定位置 |
| 给普通 gzip 文件建 GZI 失败 | 普通 gzip 没有块结构；先 `ow convert` 成 bgzip FASTA |
| 重新压缩后查询出错 | 即使碱基相同，块边界也可能改变；使用与新文件配套的 `.fai` 和 `.gzi` |
| `detect` 返回 `Bed` 或不能识别 `.gzi` | GZI 没有专门的探测分支；使用 `Gzi(path)` 显式声明 |

## 相关格式

- [FASTA](fasta.md)：bgzip FASTA 是 GZI 服务的数据文件。
- [FAI](fai.md)：与 GZI 配合定位碱基。
- [BAI](bai.md)、[CSI](csi.md)：BAM / BCF 的区域索引，使用 BGZF 虚拟偏移，不依赖 GZI。

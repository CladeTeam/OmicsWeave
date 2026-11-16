# FAI：FASTA 序列索引

[格式总览](README.md) · [FASTA](fasta.md) · [GZI](gzi.md) · [查询](../query.md)

`.fai` 是 samtools `faidx` 定义的 FASTA 索引：每条序列一行，记录序列名、长度、首个碱基在文件中的位置和换行布局。有了它，读取任意一段碱基只需一次算术和一次读取，无需扫描整个文件。

在 OmicsWeave 中，普通 FASTA 和 bgzip FASTA 都依赖 `.fai` 才能交给 `SequenceReader` 查询；bgzip FASTA 还需要 [GZI](gzi.md)。`.fai` 同时是一份染色体长度表，可直接作为转换的 `chrom_sizes`。

## 速览

| 项目 | 说明 |
|---|---|
| 作用 | 定位 FASTA 中每条序列、每个碱基的位置 |
| 命名 | `<fasta>.fai`，例如 `genome.fa.fai`、`genome.fa.gz.fai` |
| 存储 | 制表符分隔的纯文本，每条序列一行 |
| 配套文件 | 纯文本 FASTA 单独使用；bgzip FASTA 还需要 `.gzi` |
| 生成工具 | `samtools faidx`、`Fasta.build_index()`、`FastaWriter`、FASTA → bgzip 转换 |
| OmicsWeave | 建索引、查询时读取、作为 `chrom_sizes` |

## 文件结构

```text
chr1	12	6	4	5
chr2	8	27	4	5
```

| 列 | 名称 | 含义 |
|---|---|---|
| 1 | NAME | 序列名（header 中 `>` 后的第一个词） |
| 2 | LENGTH | 碱基数 |
| 3 | OFFSET | 该序列第一个碱基在文本中的字节偏移 |
| 4 | LINEBASES | 每个完整行的碱基数 |
| 5 | LINEWIDTH | 每个完整行的字节数，包含换行符（`\n` 为 1 字节，`\r\n` 为 2 字节） |

上面的索引对应下面这个每行 4 个碱基的 FASTA：

```text
>chr1          ← 6 字节，chr1 的首个碱基位于偏移 6
ACGT
Nacg
tRYK
>chr2          ← 位于 6 + 3 × 5 = 21，再加 6 字节 header，chr2 的首个碱基位于 27
TTGG
CCAA
```

## 原理：如何定位一个碱基

第 `pos` 个碱基（0-based）的字节偏移为：

```text
OFFSET + (pos // LINEBASES) × LINEWIDTH + pos % LINEBASES
```

例如 chr1 的第 5 个碱基（`a`）：`6 + (5 // 4) × 5 + 5 % 4 = 12`。读取 `[start, end)` 时，只需算出起止字节，读出这段文本并去掉其中的换行符。

这套算术成立的前提是**同一条序列除最后一行外各行长度相同**。因此 OmicsWeave 与 samtools 一样，遇到行长不一致的序列会拒绝建索引；重新排版可用 [FastaWriter](fasta.md#从内存写入) 或 `ow convert`。

对 bgzip FASTA，OFFSET 指的是**解压后**文本中的位置。FAI 给出要读的文本范围，再由 [GZI](gzi.md) 找出覆盖这段范围的压缩块。

## 建索引的规则

OmicsWeave 生成的 `.fai` 与 `samtools faidx` 逐字节一致，规则相同：

- 序列名从 `>` 后的第一个非空白字节开始，到下一个空白字符结束，其余描述被忽略；名称必须是 UTF-8。
- 一条序列的各行长度必须一致，只有最后一行可以更短，且较短的行之后只能是下一个 header 或空行。
- 碱基按可打印的非空白字节计数；`\r\n` 换行的文件也可以建索引。
- 同名序列只保留第一条；没有任何序列行的记录被丢弃（但文件最后一条记录为空时报错）。
- 普通 gzip 无法建索引（与 htslib 相同）；先转换为 bgzip。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 建索引 | `Fasta.build_index()` → `Fai`；`BgzipFasta.build_index()` → `(Fai, Gzi)`，一次扫描完成 |
| 随写随建 | `FastaWriter`、FASTA → bgzip 转换直接根据写出的布局生成，不回读文件 |
| 查询 | `SequenceReader(path, fai=...)`；默认使用 `<fasta>.fai` |
| 路径对象 | `Fai.beside(path)` 返回 `<path>.fai` 的格式对象，不创建、不校验 |
| 染色体长度 | 转换的 `chrom_sizes` / `--chrom-sizes` 接受 `.fai` |
| 保存方式 | `build_index()` 先写临时文件再重命名，中途失败不会留下半个索引 |

限制：

- reader 打开时不会自动建索引，缺少索引报 `IoError`。
- `.fai` 没有固定的文件签名，`formats.detect` 不能识别它；需要时用 `Fai("path")` 声明。

## 使用

### 创建

```python
from omicsweave.formats import BgzipFasta, Fasta

fai = Fasta("genome.fa").build_index()             # genome.fa.fai
fai, gzi = BgzipFasta("genome.fa.gz").build_index()  # genome.fa.gz.fai 和 .gzi
print(fai.path)
```

按文件类型选择其中一条。`FastaWriter` 和 `ow convert genome.fa genome.fa.gz` 会一并生成索引，无需另外调用。

### 查询时指定

```python
from omicsweave.query import SequenceReader

with SequenceReader("genome.fa", fai="index/genome.fai") as genome:
    print(genome.fetch(("chr1", 0, 100)).text)
```

未显式提供时，reader 查找 `<fasta>.fai`。`Fai("index/genome.fai")` 可声明一个已有索引的格式对象。

### 作为染色体长度表

```sh
ow convert peaks.bed peaks.bb --chrom-sizes genome.fa.fai
```

转换只使用前两列（名称和长度）。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 打开 FASTA 提示缺少 `.fai` | 调用 `build_index()`，或用 `fai=` 指定位置 |
| 建索引报行长不一致 | 重新排版 FASTA 后再建索引 |
| 查询结果错位 | FASTA 被修改或替换过，旧索引失效；重建 |
| 普通 gzip FASTA 无法建索引 | 先 `ow convert` 成 bgzip FASTA |
| `detect` 认不出 `.fai` | 正常现象：FAI 没有签名 |

## 相关格式

- [FASTA](fasta.md)：FAI 服务的数据文件。
- [GZI](gzi.md)：bgzip FASTA 需要与 FAI 配合使用的块偏移索引。
- [2bit](2bit.md)：自带目录，不需要 FAI。

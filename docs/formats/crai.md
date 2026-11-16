# CRAI：CRAM 区域索引

[格式总览](README.md) · [CRAM](cram.md) · [BAI](bai.md) · [查询](../query.md)

`.crai` 是 [CRAM](cram.md) 的**区域索引**：一个 gzip 压缩的文本文件，每行描述一个 slice 覆盖的参考范围及其在 CRAM 文件中的位置。OmicsWeave 查询时，先根据参考编号与起点在 CRAI 中找到一个起始 container，再从该 container 顺序扫描，按 container / slice 的范围跳过不重叠的数据，直到超过查询区域。它不会把所有重叠 slice 的偏移先收集成列表逐个 seek。与 BAM 的 [BAI](bai.md) 相比，它不需要 bin 层级，因为 CRAM 本身已经按 slice 组织好了数据。

在 OmicsWeave 中，`AlignmentReader` 打开 CRAM 时必须有 CRAI；可以用 `Cram.build_index()` 生成，生成过程不需要参考序列。

## 速览

| 项目 | 说明 |
|---|---|
| 作用 | 坐标排序 CRAM 的区域定位 |
| 常见后缀 | `<cram>.crai`（如 `reads.cram.crai`） |
| 存储 | gzip 压缩的 tab 分隔文本，每行 6 个整数 |
| 适用范围 | OmicsWeave 支持坐标排序的 CRAM 3.0 / 3.1 |
| 典型来源 | `samtools index reads.cram`、`Cram.build_index()` |
| OmicsWeave | 生成、自动查找、显式指定、格式识别 |

## 文件结构

解压后每行对应一个 slice（多参考 slice 中的每条参考各一行），六列依次为：

| 列 | 含义 |
|---|---|
| 1 | 参考序列编号（按 header 中 `@SQ` 的顺序，从 0 开始）；`-1` 表示未定位读段 |
| 2 | slice 中比对的起始位置，**1-based** |
| 3 | 跨度（覆盖的参考碱基数） |
| 4 | slice 所在 container 在文件中的字节偏移 |
| 5 | slice 相对 container header 末尾的偏移 |
| 6 | slice 的字节数 |

例如：

```text
0	101	4900	399	234	424
1	10	5	1076	201	86
-1	0	1	1381	179	90
```

第一行表示：第 0 条参考（如 `chr1`）的 101–5000 位置上的记录，位于文件偏移 399 处 container 中、header 之后第 234 字节开始的 424 字节 slice。最后一行是未定位读段所在的 slice。

可以直接查看：

```sh
zcat reads.cram.crai | head
```

## 关键概念

### 为什么不需要参考

CRAI 只记录 slice 的位置和覆盖范围，这些信息都在 container 和 slice header 中，不涉及读段碱基。单参考 slice 的范围直接取自 slice header；只有多参考 slice 需要解码记录的参考编号和位置，才能为每条参考分别写一行。两者都不需要参考序列。

### CRAI 不能代替参考

CRAI 帮助定位，但不含任何参考碱基。读取 CRAM 中的序列仍需要压缩时使用的参考，见 [CRAM](cram.md#参考序列依赖)。

### 与数据绑定

CRAI 中的偏移是 CRAM 文件中的字节位置。CRAM 被重写（重新排序、重新压缩、换版本）后必须重建索引。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 生成 | `Cram.build_index()`，不需要参考；输出 gzip 压缩的索引行，比较内容时应先解压，不能假定压缩字节与 `samtools index` 相同 |
| 读取与查询 | 与 htslib 的 `cram_index_load` / `cram_index_query` 相同的加载和查找方式 |
| 自动查找 | 未指定 `index` 时依次查找 `<cram>.crai`、`<去扩展名>.crai` |
| 显式指定 | `AlignmentReader(path, index="index/reads.crai")` |
| 格式识别 | `formats.detect` 把首行为 6 个 tab 分隔整数的 gzip 文本识别为 `Crai` |
| 路径对象 | `Crai.beside(path)` 构造 `<path>.crai`，不读写文件 |

缺少 CRAI 时 reader 报 `IoError`，不会自动创建。

## 创建

```python
from omicsweave.formats import Cram

index = Cram("reads.cram").build_index()
print(index.path)  # reads.cram.crai
```

输入需要按坐标排序。

## 使用

下例中的 `genome.fa.gz` 必须是 BGZF FASTA，且配有 `.fai` / `.gzi`，普通 gzip FASTA 不能作区域读取参考。可先准备：

```sh
bgzip -c genome.fa > genome.fa.gz
```

```python
from omicsweave.formats import BgzipFasta

BgzipFasta("genome.fa.gz").build_index()
```

```python
from omicsweave.query import AlignmentReader

with AlignmentReader(
    "reads.cram", index="index/reads.crai", reference="genome.fa.gz"
) as reads:
    print(reads.index)  # Crai('index/reads.crai')
    hits = reads.fetch(("chr1", 100, 200), "all")
```

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 提示缺少索引 | `Cram(path).build_index()`，或把索引放在 `<cram>.crai` |
| 有 CRAI 仍无法读取序列 | CRAI 不含参考；提供压缩时的参考 |
| 能否用 BAI / CSI | 不能，CRAM 只使用 CRAI |
| CRAM 改写后结果异常 | 旧 CRAI 失效，重建索引 |

## 相关格式

- [CRAM](cram.md)：CRAI 所索引的数据文件。
- [BAI](bai.md)、[CSI](csi.md)：BAM / BCF 的 binning 索引。

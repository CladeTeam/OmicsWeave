# CSI：BAM / BCF 区域索引

[格式总览](README.md) · [BCF](bcf.md) · [BAM](bam.md) · [BAI](bai.md) · [查询](../query.md)

CSI（Coordinate-Sorted Index）是坐标排序的 BGZF 文件的**区域索引**：给定一个基因组区间，它告诉读取方应该从压缩文件的哪些位置开始解压，才能拿到所有与该区间重叠的记录。它本身不含任何比对或变异数据，必须与建索引时的那个数据文件一起使用。

在 OmicsWeave 中，CSI 是 [BCF](bcf.md) 唯一的区域索引，也是 [BAM](bam.md) 的两种区域索引之一（另一种是 [BAI](bai.md)）。`VariantReader`、`ConsensusReader` 和 `AlignmentReader` 都通过它按区域读取；`Bam.build_index`、`Bcf.build_index` 和转换时的 `index=` 都可以生成它。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 区域 → 压缩文件中的虚拟偏移区间 |
| 常见后缀 | `.csi`，通常命名为 `<数据文件>.csi` |
| 存储 | BGZF 压缩的二进制，解压后以 magic `CSI\1` 开头 |
| 适用数据 | 按坐标排序、BGZF 压缩的 BCF、BAM（以及 tabix 索引的文本文件） |
| 坐标范围 | 由 `min_shift` 和层数决定，可覆盖远超 512 Mb 的染色体 |
| 规范 | hts-specs CSIv1 |
| 典型来源 | `bcftools index`、`samtools index -c`、`tabix -C`、`ow convert ... --index` |
| OmicsWeave | 为 BAM / BCF 生成，供 BAM、BCF 区域查询和 consensus 使用 |

## 格式由来

htslib 最初为 BAM 设计了 BAI 索引，它采用固定的分箱方案：最小 bin 16 kb、共 6 层，因此单条参考序列最长只能到 2²⁹ = 536,870,912 bp。小麦、大麦、蝾螈等物种的染色体超过这个长度，BAI 就无法使用。

CSI 由 htslib 作者在 BAI 的基础上推广而来：分箱方式相同，但最小 bin 大小（`min_shift`）和层数（`depth`）可以配置，并写在索引文件中。这样既能索引超长染色体，也能按需调整索引粒度。BCF 从一开始就使用 CSI。

## 索引如何工作

### 虚拟偏移

BGZF 文件由许多独立压缩的小块组成。文件中一条记录的位置用 64 位**虚拟偏移**表示：

```text
虚拟偏移 = (压缩块在文件中的字节位置 << 16) | 记录在解压后块内的偏移
```

高 48 位用于找到块并单独解压它，低 16 位定位块内的记录。索引中保存的全部是这种偏移。

### 分层分箱

CSI 把每条参考序列切成层级嵌套的 bin：

```text
层 0：            [──────────────── 一个 bin 覆盖 2^(min_shift + 3·depth) ────────────────]
层 1：            [────────] [────────] [────────] ... 共 8 个
层 2：            每个父 bin 再分为 8 个子 bin
...
层 depth：        最小的 bin，每个覆盖 2^min_shift 个碱基
```

每条记录放入**能完整容纳它**的最小 bin：短记录落在底层的小 bin，跨越边界的长记录（长删除、长读段、结构变异）落在上层的大 bin。bin 按层连续编号，第 `l` 层的第一个编号是 `(8^l − 1) / 7`。

每个 bin 记录一组 **chunk**：该 bin 中的记录在文件中占据的虚拟偏移区间 `[起, 止)`。CSI 还为每个 bin 保存一个 `loffset`，表示该 bin 范围内第一条记录可能出现的最小偏移，用于跳过不可能重叠的 chunk。BAI 用单独的“线性索引”做同样的事。

每条参考序列还有一个**伪 bin**，记录已比对和未比对记录的计数；文件末尾记录没有坐标的记录数。

### 一次区域查询

1. 计算与区间 `[beg, end)` 重叠的所有 bin（每层若干个）。
2. 收集这些 bin 的 chunk，去掉止点早于 `loffset` 下界的部分，排序并合并相邻 chunk。
3. 依次 seek 到每个 chunk，解压并解析记录，丢弃不与区间重叠的那些。

所以索引只能缩小读取范围，最终是否重叠仍由记录本身的坐标判断。

### min_shift 与 depth

| 参数 | 含义 | 默认 |
|---|---|---|
| `min_shift` | 最小 bin 为 `2^min_shift` bp | 14（16 kb），与 BAI 相同 |
| `depth` | 底层之上的层数 | 按最长参考序列自动确定 |

`min_shift` 越小，底层 bin 越细，小区域查询读得越少，但索引越大；一般保持默认 14 即可。`depth` 需满足 `2^(min_shift + 3·depth)` 不小于最长的参考序列。BAI 相当于固定 `min_shift=14, depth=5` 的 CSI。

### 与 BAI、tabix 的区别

| | BAI | CSI | tabix `.tbi` |
|---|---|---|---|
| 适用 | BAM | BAM、BCF、tabix 文本 | bgzip 文本（VCF、BED、GFF 等） |
| 分箱 | 固定 14 / 5 层 | 可配置 | 固定，同 BAI |
| 最长参考 | 2²⁹ bp | 由参数决定 | 2²⁹ bp |
| 跳过 chunk | 线性索引 | 每个 bin 的 `loffset` | 线性索引 |
| 文件压缩 | 通常不压缩 | BGZF | BGZF |

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 格式识别 | 按解压后的 magic `CSI\1` 识别为 `Csi` |
| 为 BCF 建索引 | `Bcf.build_index(min_shift=14, threads=None)`，等同 `bcftools index` |
| 为 BAM 建索引 | `Bam.build_index(kind="csi", min_shift=14)`；`kind="auto"` 在参考超过 2²⁹ bp 时自动选择 CSI |
| 转换时建索引 | VCF → BCF：`index=True`（`min_shift` 固定为 14）；SAM → BAM、BAM 排序、CRAM → BAM：`index="csi"`，可用 `min_shift` |
| 层数 | 按最长参考序列自动计算；BCF header 中没有任何 contig 长度时按 2³¹ − 1 计算 |
| 读取 | `AlignmentReader`（BAM）、`VariantReader`、`ConsensusReader` |
| 查找顺序 | BAM：`<bam>.csi`、去扩展名的 `.csi`、`<bam>.bai`、去扩展名的 `.bai`；BCF：`<bcf>.csi`、去扩展名的 `.csi` |
| 自定义位置 | 各 reader 的 `index=` 参数 |
| 查询语义 | 与 htslib `hts_itr_query` 一致，返回的记录集合与 `sam_itr_queryi` / `bcf_itr_queryi` 相同 |

限制：

- 打开 reader 时不会自动建索引，缺少索引时报错并提示建索引的方法。
- 建索引要求记录已按坐标排序；CSI 不会替数据排序。
- 不读取 tabix 格式的文本索引；文本 VCF 需先转为 BCF。
- `min_shift` 在转换中只对 BAM 的 CSI 生效；VCF → BCF 的索引使用默认 14，需要其他值时转换后再调用 `Bcf.build_index(min_shift=...)`。

## 创建

```python
from omicsweave.formats import Bam, Bcf

Bcf("cohort.bcf").build_index()                     # cohort.bcf.csi
Bam("reads.bam").build_index(kind="csi")            # reads.bam.csi
Bam("reads.bam").build_index(kind="csi", min_shift=12, threads=4)
Bam("wheat.bam").build_index()                      # auto：参考超过 2^29 bp 时写 CSI
```

`threads` 默认使用进程可用的全部核心解压数据文件。对超长参考强制 `kind="bai"` 会报错：`a reference of ... bases does not fit a BAI index; use CSI`。

在转换的同时生成：

```sh
ow convert cohort.vcf.gz cohort.bcf --index
ow convert unsorted.bam sorted.bam --sort --index csi
ow convert unsorted.bam sorted.bam --sort --index csi --min-shift 12
ow convert reads.cram reads.bam --reference genome.fa.gz --index csi
```

```python
from omicsweave import convert

report = convert.run("unsorted.bam", "sorted.bam", sort=True, index="csi")
print(report.companions)  # (Csi('sorted.bam.csi'),)
```

## 使用

```python
from omicsweave.query import AlignmentReader, VariantReader

with AlignmentReader("reads.bam") as reads:
    print(reads.index)                   # 实际使用的索引，如 Csi('reads.bam.csi')

with VariantReader("cohort.bcf", index="index/cohort.csi") as variants:
    gt = variants.genotypes(("chr1", 0, 100_000))
```

按常规命名放在数据文件旁边时会被自动找到，其他位置通过 `index=` 传入。reader 的 `index` 属性给出实际使用的索引文件，可用于确认读到的是 CSI 还是 BAI。

`Csi.beside(path)` 只构造 `<path>.csi` 的格式对象，不创建或校验文件。

## 注意事项

- **索引与数据绑定**：索引保存的是数据文件中的字节位置。重新压缩、排序、改写数据文件后，即使内容相同也必须重建索引；把旧索引与新文件混用会读到错误的记录或报错。
- **排序要求**：记录必须按参考序列在 header 中的顺序、再按起点排序。VCF → BCF、SAM → BAM、CRAM → BAM 保留输入顺序，未排序输入不能同时建索引。
- **quant**：`quant` 顺序扫描主输入 BAM，不需要它的索引；但 `Task.variants`、`Alleles` 等 BCF 辅助输入需要配套 CSI。
- **其他格式**：CRAM 使用 [CRAI](crai.md)，FASTA 使用 [FAI](fai.md) + [GZI](gzi.md)，bigWig / bigBed / `.hic` 自带内置索引，都不使用 CSI。quant 为 `fragments.tsv.gz` 等表格输出的 `.tbi` 是 tabix 索引，供下游工具使用。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| `no index of ... found` | 没有索引；调用 `build_index()`，或用 `index=` 指定位置 |
| 建索引报 `unsorted positions` | 数据未排序；BAM 用 `ow convert x.bam sorted.bam --sort`，VCF 先排序再转换 |
| BAI 报参考过长 | 染色体超过 2²⁹ bp；改用 `kind="csi"` |
| `option min_shift does not apply to VCF to BCF` | VCF → BCF 不接受 `min_shift`；转换后用 `Bcf.build_index(min_shift=...)` |
| 同时有 `.csi` 和 `.bai` | BAM reader 优先使用 `.csi`；需要另一个时用 `index=` 指定 |

## 相关格式

- [BAI](bai.md)：BAM 的固定分箱索引，CSI 的前身。
- [BCF](bcf.md)、[BAM](bam.md)：使用 CSI 的数据格式。
- [CRAI](crai.md)：CRAM 的区域索引。

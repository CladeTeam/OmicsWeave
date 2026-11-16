# BAI：BAM 区域索引

[格式总览](README.md) · [BAM](bam.md) · [CSI](csi.md) · [查询](../query.md)

`.bai` 是坐标排序 [BAM](bam.md) 的**区域索引**：它把基因组划分成层级化的 bin，记录每个 bin 中的记录在 BAM 文件里从哪里开始、到哪里结束。有了它，查询 `chr1:100-200` 时只需读取几个 BGZF 块，而不必从头扫描整个 BAM。BAI 只包含定位信息，不包含比对本身，必须与对应的 BAM 一起使用。

在 OmicsWeave 中，`AlignmentReader` 通过 BAI 或 [CSI](csi.md) 读取 BAM 区域；索引可以用 `Bam.build_index()` 生成，也可以在 SAM → BAM、BAM 排序时一并写出。[quant](../quant.md) 顺序扫描 BAM，不需要索引。

## 速览

| 项目 | 说明 |
|---|---|
| 作用 | 坐标排序 BAM 的区域定位 |
| 常见后缀 | `<bam>.bai`（如 `reads.bam.bai`），也常见 `reads.bai` |
| 存储 | 二进制，未压缩 |
| 适用范围 | 每条参考序列不超过 `2**29` bp（约 512 Mb） |
| 典型来源 | `samtools index`、`ow convert --index`、`Bam.build_index()` |
| OmicsWeave | 生成、自动查找、显式指定 |

## 格式由来

BAI 随 BAM 一起在 samtools 中引入，采用 UCSC 基因组浏览器的 **binning 方案**：把每条参考序列按固定的层级切成 bin，每条记录放进能完整包含它的最小 bin。查询时只需检查与区域重叠的少量 bin。之后提出的 [CSI](csi.md) 是它的推广，取消了长度限制。

## 文件结构

```text
magic "BAI\1"
n_ref
对每条参考序列：
    n_bin
    对每个 bin：bin 编号、n_chunk、chunk 列表（起止虚拟偏移）
    n_intv
    线性索引：每 16 kb 窗口一个虚拟偏移
（可选）n_no_coor：没有坐标的读段数
```

### bin 层级

BAI 固定使用 6 层 bin，最细一层是 16 kb（`2**14`），每上一层放大 8 倍：

| 层 | bin 大小 | bin 数 | 编号范围 |
|---|---|---|---|
| 0 | 512 Mb（`2**29`） | 1 | 0 |
| 1 | 64 Mb | 8 | 1–8 |
| 2 | 8 Mb | 64 | 9–72 |
| 3 | 1 Mb | 512 | 73–584 |
| 4 | 128 kb | 4096 | 585–4680 |
| 5 | 16 kb | 32768 | 4681–37448 |

一条记录放入**能完整包含其参考跨度的最小 bin**。短读段通常落在第 5 层；跨越 16 kb 边界的读段、长读长或跨长内含子的 RNA-seq 读段会落到更上层。BAM 记录中的 `bin` 字段保存的就是这个编号。

这是记录的初始分箱规则。写出索引前，OmicsWeave 会把压缩文件跨度较小、且父 bin 已存在的子 bin 的 chunks 并入父 bin，再合并相邻 chunks。因此索引中最终保存的 bin 不一定与每条 BAM 记录的 `bin` 字段一一对应。

### chunk 与虚拟偏移

每个 bin 记录一组 **chunk**：一个 chunk 是 BAM 文件中一段连续的记录，用起止[虚拟偏移](bam.md#bgzf-分块压缩)表示。查询一个区域时：

1. 计算与区域重叠的所有 bin（每层若干个）。
2. 收集这些 bin 的 chunk。
3. 用线性索引排除一定在区域之前结束的 chunk。
4. 读取剩下的 chunk，逐条检查记录是否真正与区域重叠。

### 线性索引

线性索引为每个 16 kb 窗口记录“与该窗口重叠的记录中最小的虚拟偏移”。上层 bin 可能很大（例如 bin 0 覆盖整条染色体），线性索引能让查询跳过这些 bin 中远在区域之前的记录。

### 元数据伪 bin

有记录的参考序列还有一个编号为 37450 的伪 bin，记录该参考序列的数据范围以及已比对 / 未比对的读段数。OmicsWeave 为没有记录的参考写空索引项，不写伪 bin。`samtools idxstats` 使用这类元数据统计记录数。

## 关键概念

### 2^29 的限制

bin 0 覆盖 `2**29` = 536,870,912 bp，所以 BAI 无法索引更长的参考序列。人和小鼠的染色体都在此范围内，但一些植物和两栖动物基因组（例如小麦、蝾螈的部分染色体）超出了，需要使用 [CSI](csi.md)。

OmicsWeave 的 `kind="auto"`（默认）会检查 header 中所有参考序列的长度：都不超过 `2**29` 时写 BAI，否则写 CSI。

这个上限对应 0-based `[0, 2**29)`，因此区间右端恰好为 `2**29` 可以索引。tabix 的 `.tbi` 使用相同固定分箱上限，但它是文本索引，不能与 `.bai` 互换。quant 的 fragments / cytosines 表在所有参考长度均不超过此上限时写 `.tbi`，否则写带 BED 配置的 [CSI](csi.md#quant-的-tabix-csi)。

### 索引与数据绑定

BAI 的 chunks 和线性索引保存的是 BGZF 虚拟偏移：高 48 位为压缩块在 BAM 文件中的字节偏移，低 16 位为解压后块内偏移。BAM 一旦重新排序、重新压缩或以任何方式改写，旧索引就失效，必须重建；即使内容相同，压缩块的位置也可能变了。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 生成 | `Bam.build_index("bai")` 或 `kind="auto"`；转换时 `--index` / `index=True` |
| 内容 | 与 `samtools index` 相同的 bins、chunks 和线性偏移 |
| 自动查找 | 未指定 `index` 时依次查找 `<bam>.csi`、`<去扩展名>.csi`、`<bam>.bai`、`<去扩展名>.bai` |
| 显式指定 | `AlignmentReader(path, index="index/reads.bai")` |
| 格式识别 | `formats.detect` 按文件开头的 `BAI\1` magic 识别为 `Bai`；gzip 压缩的 BAI 不走此分支，会回退为 `GzipBed` |
| 路径对象 | `Bai.beside(path)` 构造 `<path>.bai`，不读写文件 |

打开 reader 不会自动生成索引；缺少索引时报 `IoError`。

## 创建

```python
from omicsweave.formats import Bam

index = Bam("reads.bam").build_index(kind="bai")
print(index.path)  # reads.bam.bai
```

输入必须已按坐标排序，否则报错。`threads` 控制解压 BAM 的线程数，默认使用进程可用的全部核心。

在转换时一并生成：

```sh
ow convert sorted.sam sorted.bam --index bai
ow convert unsorted.bam sorted.bam --sort --index
```

## 使用

```python
from omicsweave.query import AlignmentReader

with AlignmentReader("reads.bam", index="index/reads.bai") as reads:
    print(reads.index)  # Bai('index/reads.bai')
    hits = reads.fetch(("chr1", 100, 200))
```

自动查找时 CSI 优先于 BAI，这与 htslib 一致；同时存在两种索引时，可以用 `index=` 明确选择。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 提示缺少索引 | `Bam(path).build_index()`，或把索引放在 `<bam>.bai` |
| 参考序列太长无法写 BAI | 使用 CSI：`build_index("csi")` 或 `--index csi` |
| 排序后查询结果异常 | 旧索引失效，重建索引 |
| CRAM / BCF 能否用 BAI | 不能：CRAM 用 [CRAI](crai.md)，BCF 用 [CSI](csi.md) |

## 相关格式

- [BAM](bam.md)：BAI 所索引的数据文件，虚拟偏移的定义。
- [CSI](csi.md)：BAI 的推广，支持更长的参考序列，也用于 BCF。
- [CRAI](crai.md)：CRAM 的索引。

# bedGraph

[格式总览](README.md) · [查询](../query.md) · [转换](../convert.md) · [bigWig](bigwig.md) · [wig](wig.md) · [BED](bed.md)

bedGraph 用四列文本表示**基因组上的数值信号**：染色体、起点、右端和这段区间上的值。它是覆盖度、富集度、甲基化比例等轨道最常见的文本中间格式，许多工具（`bedtools genomecov -bg`、MACS2 的 `--bdg` 输出、deepTools `--outFileFormat bedgraph`）都直接产出 bedGraph。

在 OmicsWeave 中，bedGraph 是**转换源**：转为 [bigWig](bigwig.md) 后用 `TrackReader` 读取逐碱基信号和分箱统计；也可以转为 [bigBed](bigbed.md)，把每行当作带值的区间记录保存。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 区间 → 数值，每个位置最多一个值 |
| 常见后缀 | `.bedGraph`、`.bedgraph`、`.bdg`、`.bg`，压缩后加 `.gz` |
| 存储 | 纯文本，可用 gzip / bgzip 压缩 |
| 坐标 | 0-based 半开，与 BED 相同 |
| 随机访问 | 不支持；转为 bigWig 后支持 |
| 典型来源 | `bedtools genomecov -bg`、MACS2 / MACS3、deepTools、UCSC 下载 |
| OmicsWeave | 格式识别、转换为 bigWig 或 bigBed |

## 格式由来

bedGraph 由 UCSC 基因组浏览器定义，用来上传连续数值轨道。它借用了 BED 的前三列和坐标约定，再加一列数值，因此既能被 BED 工具链处理，又能表达信号。UCSC 的 `bedGraphToBigWig` 把它压缩成 bigWig；OmicsWeave 的转换按 `bedGraphToBigWig` 的方式读取文本。

和 [wig](wig.md) 相比，bedGraph 每行都写出完整坐标，因此可以表达长度不一的区间，适合稀疏或分段常数的信号；wig 的 `fixedStep` 在信号密集、步长固定时更紧凑。

## 文件结构

```text
track type=bedGraph name="ATAC" description="ATAC-seq coverage"
chr1	100	150	2.5
chr1	150	160	4
chr1	200	300	8.0
chr2	0	1000	0.25
```

| 列 | 含义 |
|---|---|
| 1 | 染色体名 |
| 2 | 起点，0-based，包含 |
| 3 | 右端，不包含 |
| 4 | 数值，整数或浮点 |

可选的 header 行：

- `track type=bedGraph ...`：UCSC 轨道定义，`type=bedGraph` 声明这是 bedGraph。
- `browser ...`：浏览器显示设置。
- `#` 开头的注释。

上例中 `chr1 100 150 2.5` 表示第 100–149 号碱基（0-based）的值都是 2.5。`chr1 160–200` 没有任何一行覆盖，表示**没有数据**。

## 关键概念

### 分段常数

bedGraph 描述的是“分段常数”函数：一行内所有碱基取同一个值。相邻且值相同的区间可以合并成一行，结果不变；反过来，一个长区间也可以拆成多行。

### 没有数据 ≠ 0

未被任何行覆盖的位置是“没有数据”。很多工具输出覆盖度时会省略 0 值区间，有的则会显式写出 0（如 `bedtools genomecov -bga`）。转为 bigWig 后，前者在 `TrackReader.signal` 中是 NaN，后者是 0。是否把 NaN 当作 0 需要依据数据来源判断。

### 不能重叠

信号在一个位置上只能有一个值，所以区间不能重叠。两个工具输出合并、或不同链的信号写进同一文件时容易产生重叠，需要先拆分或合并。

### 与 BED4 的区别

bedGraph 和 BED4 都是四列，区别在第四列的含义：bedGraph 是数值，BED4 是名称。数字形式的“名称”（如 peak 编号）会让 BED4 看起来像 bedGraph，见下文的识别规则。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 文件对象 | `BedGraph`（纯文本）、`GzipBedGraph`（gzip 或 bgzip） |
| 格式识别 | 开头最多 64 KiB 中有 `track type=bedGraph`，或样本数据行都是四列、坐标可解析为 `uint32`、末列可解析为浮点数 |
| 区域查询 | 不直接支持；`query.open` 会报错并给出转换命令 |
| 转为 bigWig | 存储数值信号，支持逐碱基查询和分箱统计 |
| 转为 bigBed | 保留四列区间语义，数值作为附加字段文本 |

## 转为 bigWig

```sh
ow convert signal.bedGraph signal.bw --chrom-sizes genome.chrom.sizes
ow convert signal.bedGraph.gz signal.bw --chrom-sizes genome.fa.fai -t 4 -m 2G
```

```python
from omicsweave import convert

report = convert.run("signal.bedGraph", "signal.bw", chrom_sizes="genome.2bit")
print(report.records)  # 写入的信号区间数
```

转换规则：

- **输入**：纯文本、gzip、BGZF；`track`、`browser`、`#` 开头的行和空行跳过。字段以空白分隔。
- **排序**：不要求有序，按染色体（首次出现的顺序）和起点排序；超出内存预算时外排序。
- **区间**：必须非空（`start < end`），排序后相邻区间不能重叠，重叠时报错并指出两行的行号。
- **数值**：按双精度解析，再转为 `float32`；转换后必须有限，`nan`、`inf` 或溢出 float32 的值报错。小数可能舍入或下溢，例如 `1e-300` 会存为 0；不做额外缩放或平滑。
- **染色体长度**：`chrom_sizes` 接受两列长度表、`.fai` 或 `.2bit`。提供时，不在长度表中的染色体和超出长度的区间都会报错，不会自动裁剪；省略时用每条染色体最大的右端作为长度。
- 输出 bigWig 默认尝试生成 zoom 层，只保留比上一层记录更少的层，因此小文件可能没有 zoom 层。块与 zoom 布局可用 `bbi=BbiOptions(...)` 或 `--block-size` / `--items-per-slot` / `--zoom-levels` 调整；`--sidecar` 可同时写出 [sidecar](sidecar.md)。

转换后读取：

```python
from omicsweave.query import TrackReader

with TrackReader("signal.bw") as track:
    print(track.signal(("chr1", 100, 300)).values)                # 逐碱基，空隙为 NaN
    print(track.bins(("chr1", 100, 300), 100, "mean")["mean"])    # 分箱均值
    print(track.intervals(("chr1", 100, 300)).value)              # 存储的原始区间
```

`intervals` 返回的正是 bedGraph 中的行（裁剪到查询区域），可用于检查转换结果。

## 转为 bigBed

需要把每行当作“带值的区间记录”保存，例如保留原始文本或与其他区间一起做重叠查询时，可以转为 bigBed：

```sh
ow convert signal.bedGraph signal.bb --chrom-sizes genome.chrom.sizes
```

```python
from omicsweave.query import IntervalReader

with IntervalReader("signal.bb") as records:
    hits = records.fetch(("chr1", 0, 1000))
    print(hits.start, hits.end, hits.annotation.to_list())  # 数值以文本保存，如 ['2.5', '8.0']
```

这条路由按 [BED](bed.md) 的规则处理（允许重叠，数值不做浮点校验）。需要信号数组或分箱统计时，目标应选择 bigWig。

## 识别与显式指定

`formats.detect` 读取前 64 KiB 判断格式。以下情况会被识别为 bedGraph：

- `track` 行写有 `type=bedGraph`；
- 样本中有数据行，每个完整数据行都恰好四列，起点和右端可解析为 `uint32`，最后一列可解析为浮点数。

如果一个四列文件实际上是 BED4（第四列是数字编号的名称），显式按 BED 转换：

```python
from omicsweave.formats import Bed

convert.run(Bed("regions.txt"), "regions.bb", chrom_sizes="genome.chrom.sizes")
```

反过来，一个后缀为 `.txt` 的 bedGraph 也会被正确识别，文件名不参与判断。也可以用 `BedGraph(path)` 或 CLI 的 `--from bedgraph` 明确指定。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| `... overlaps ... of line N` | 区间重叠；按位置合并（求和或取均值，取决于数据含义）后再转换 |
| `empty item` | 有 `start == end` 的空区间；删除这些行 |
| `value ... is not a finite 32-bit float` | 存在 `nan`、`inf` 或超大数值；先清理 |
| `chromosome ... is not in the chromosome sizes` | 长度表与文件的染色体命名或参考版本不一致 |
| 查询结果大量 NaN | 原文件省略了无信号区间；按数据含义决定是否填 0 |
| 四列 BED 被当成 bedGraph | 用 `Bed(path)` 或 `--from bed` 显式指定 |

## 相关格式

- [bigWig](bigwig.md)：bedGraph 的带索引二进制形式，查询数值信号使用它。
- [wig](wig.md)：另一种信号文本格式，允许用固定步长紧凑表示。
- [BED](bed.md)：同样的坐标约定，表示区间而不是信号。
- [bigBed](bigbed.md)：bedGraph 也可作为区间记录转入。

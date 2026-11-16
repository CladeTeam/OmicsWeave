# wig

[格式总览](README.md) · [查询](../query.md) · [转换](../convert.md) · [bigWig](bigwig.md) · [bedGraph](bedgraph.md)

wig（wiggle）是表示**基因组数值信号**的文本格式。与每行写完整坐标的 [bedGraph](bedgraph.md) 不同，wig 先用一行声明“从哪里开始、每隔多少碱基一个值”，后面的行只写数值，因此对密集、等间距的信号（保守性分数、GC 含量、固定窗口覆盖度）非常紧凑。UCSC 的 phyloP / phastCons 保守性轨道就以 wig 形式发布。

在 OmicsWeave 中，wig 是**转换源**：转为 [bigWig](bigwig.md) 后，用 `TrackReader` 按区域读取逐碱基信号和分箱统计。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 基因组数值信号 |
| 常见后缀 | `.wig`、`.wig.gz` |
| 存储 | 纯文本，可用 gzip / bgzip 压缩 |
| 坐标 | 段声明中的 `start` 和 variableStep 位置为 **1-based** |
| 随机访问 | 不支持；转为 bigWig 后支持 |
| 典型来源 | UCSC 保守性轨道（phyloP、phastCons）、GC 含量、早期 ChIP-seq 工具 |
| OmicsWeave | 格式识别、转换为 bigWig |

## 格式由来

wig 是 UCSC 基因组浏览器最早支持的数值轨道格式之一。为了让一整条染色体的逐碱基分数尽量小，它采用“声明 + 数值”的有状态写法：声明行设定当前段的染色体、起点、步长和跨度，之后的每个数值依次落在推算出的位置上。bigWig 正是为了替代大型 wig 文件而设计的；UCSC 的 `wigToBigWig` 完成这一转换，OmicsWeave 按同样的规则读取 wig。

## 文件结构

一个 wig 文件由若干**段**组成，每段以声明行开头。常见的两种段：

### fixedStep：等间距

```text
fixedStep chrom=chr1 start=101 step=10 span=10
2.5
3.0
4.5
```

| 参数 | 含义 | 是否必需 |
|---|---|---|
| `chrom` | 染色体名 | 必需 |
| `start` | 第一个值的起点，1-based | 必需 |
| `step` | 相邻两个值的起点间隔 | 可选，默认 1 |
| `span` | 每个值覆盖的碱基数 | 可选，默认 1 |

第 i 个值（从 0 数）覆盖 1-based 闭区间 `[start + i·step, start + i·step + span - 1]`，换成 0-based 半开坐标是 `[start - 1 + i·step, start - 1 + i·step + span)`。上例三个值分别覆盖 `[100, 110)`、`[110, 120)`、`[120, 130)`。

### variableStep：显式位置

```text
variableStep chrom=chr2 span=5
201 4.0
301 7.0
```

| 参数 | 含义 | 是否必需 |
|---|---|---|
| `chrom` | 染色体名 | 必需 |
| `span` | 每个值覆盖的碱基数 | 可选，默认 1 |

数据行是“1-based 位置 + 值”。上例两行分别覆盖 0-based 的 `[200, 205)` 和 `[300, 305)`。适合位置不规则但跨度一致的信号。

### 其他行

- `track type=wiggle_0 name=... description=...`：UCSC 轨道定义，`type=wiggle_0` 声明 wig。
- `browser ...`、`#` 注释、空行。
- bedGraph 行（`chrom start end value`，0-based 半开）：一些文件在 wig 中混用 bedGraph 行，它结束当前的 step 段。

### 同一信号的三种写法

下面三段描述的是同一份信号——chr1 的 0-based `[100, 130)` 上每 10 bp 一个值：

```text
fixedStep chrom=chr1 start=101 step=10 span=10
2.5
3.0
4.5
```

```text
variableStep chrom=chr1 span=10
101 2.5
111 3.0
121 4.5
```

```text
chr1	100	110	2.5
chr1	110	120	3.0
chr1	120	130	4.5
```

## 关键概念

### 1-based 的段坐标

wig 的 `start` 和 variableStep 位置是 1-based，而 bedGraph、BED、bigWig 和 OmicsWeave 查询都是 0-based 半开。转换时 OmicsWeave 自动减一，因此 `start=101` 的第一个值在查询坐标中位于 100。混用 bedGraph 行时，这些行仍按 0-based 解释。

### step 与 span

`step` 决定值与值之间的间隔，`span` 决定每个值覆盖多宽：

- `span == step`：值首尾相接，信号连续。
- `span < step`：值之间留有空隙，空隙中没有数据。
- `span > step`：相邻值互相重叠。bigWig 每个位置只能有一个值，这种写法在转换时报错。

### 有状态的解析

数值行本身不含坐标，其位置完全取决于前面最近的声明。因此 wig 不能随意截取中间的行，也不能简单地按行排序；转换器必须从头按顺序理解段结构。OmicsWeave 在并行解析时会处理跨批次的段延续，结果与顺序读取一致。

### 没有数据 ≠ 0

段之外、或 `span < step` 留下的空隙都是没有数据，转为 bigWig 后查询得到 NaN。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 文件对象 | `Wig`（纯文本）、`GzipWig`（gzip 或 bgzip） |
| 格式识别 | `track` 行写有 `type=wiggle_0`，或出现 `fixedStep` / `variableStep` 声明 |
| 段类型 | fixedStep、variableStep，以及混用的 bedGraph 行 |
| 区域查询 | 不直接支持；`query.open` 会报错并给出转换命令 |
| 转换 | wig → bigWig |

限制：

- 不保留原始段结构：bigWig 中保存的是等价的“区间 → 值”，查询结果相同，但无法还原 fixedStep / variableStep 的文本写法。
- 当前不提供 bigWig → wig 的导出。

## 转为 bigWig

```sh
ow convert signal.wig signal.bw --chrom-sizes genome.chrom.sizes
ow convert phyloP.wig.gz phyloP.bw --chrom-sizes genome.2bit -t 8 -m 4G
```

```python
from omicsweave import convert

report = convert.run("signal.wig", "signal.bw", chrom_sizes="genome.fa.fai")
print(report.records)  # 写入的信号区间数
```

转换规则：

- **输入**：纯文本、gzip、BGZF；`track`、`browser`、`#` 开头的行和空行跳过。
- **声明**：`fixedStep` 必须有 `chrom` 和 `start`，`variableStep` 必须有 `chrom`；`step`、`span` 须为正整数，未知的设置项报错。
- **数据行**：fixedStep 段内每行一个值，variableStep 段内每行“位置 值”；列数不符合所在段时报错。
- **排序**：段可以按任意顺序出现，转换按染色体（首次出现的顺序）和起点排序。
- **重叠**：排序后相邻区间不能重叠；`span > step` 或不同段覆盖同一位置时报错。
- **数值**：以 `float32` 存储，必须是有限值；按原样保存。
- **染色体长度**：`chrom_sizes` 接受两列长度表、`.fai` 或 `.2bit`。提供时，不在表中的染色体和超出长度的区间报错；省略时以每条染色体最大的右端作为长度。
- `--sidecar` 可同时写出 [sidecar](sidecar.md)；块与 zoom 布局选项与 [bedGraph](bedgraph.md) 转换相同。

转换后读取：

```python
from omicsweave.query import TrackReader

with TrackReader("signal.bw") as track:
    values = track.signal(("chr1", 100, 130)).values   # fixedStep 示例：10 个 2.5、10 个 3.0、10 个 4.5
    runs = track.intervals(("chr1", 100, 130))
    print(runs.start, runs.end, runs.value)            # [100 110 120] [110 120 130] [2.5 3. 4.5]
```

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 报错区间重叠 | `span > step`，或多个段覆盖同一位置；调整 span 或拆分文件 |
| `fixedStep without start=` | fixedStep 声明缺少 `start` |
| `unknown setting` | 声明行有不支持的设置项（只认 `chrom`、`start`、`step`、`span`） |
| 结果整体偏移 1 bp | 自行解析 wig 时忘了 1-based；OmicsWeave 转换已自动处理 |
| 空隙中是 NaN | `span < step` 或段之间没有数据；按数据含义决定是否填 0 |

## 相关格式

- [bigWig](bigwig.md)：wig 的带索引二进制形式，查询时使用。
- [bedGraph](bedgraph.md)：每行写完整坐标的信号文本格式，可与 wig 混用。

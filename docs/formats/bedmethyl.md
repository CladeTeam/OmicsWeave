# bedMethyl

[格式总览](README.md) · [转换](../convert.md) · [定量](../quant.md) · [bigBed](bigbed.md) · [bigWig](bigwig.md) · [BED](bed.md)

bedMethyl 是记录**逐位点碱基修饰（如 DNA 甲基化）统计**的 BED 扩展格式：前 9 列是标准 BED9，后面跟着这个位点的测序覆盖度、修饰比例以及各类读段计数。它是 ENCODE 全基因组亚硫酸氢盐测序（WGBS）的发布格式，也是 Oxford Nanopore `modkit pileup` 的默认输出，因此在短读长和长读长甲基化分析中都很常见。

在 OmicsWeave 中，bedMethyl 是**转换源**，有两种去向：转为 [bigBed](bigbed.md) 保留所有列，用 `IntervalReader` 按区域查询完整记录；或挑出一个数值列（修饰比例、覆盖度、某类计数）转为 [bigWig](bigwig.md)，作为模型可直接使用的信号轨道。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 每个位点（通常是一个 C 或 A）的修饰调用统计 |
| 常见后缀 | `.bed`、`.bedmethyl`、`.bed.gz`、`.bedmethyl.gz` |
| 存储 | 制表文本，可用 gzip / bgzip 压缩 |
| 坐标 | 0-based 半开，单碱基位点为 `[pos, pos + 1)` |
| 方言 | modkit pileup（18 列）、ENCODE WGBS（11 列） |
| 随机访问 | 不支持；转为 bigBed / bigWig 后支持 |
| 典型来源 | `modkit pileup`、ENCODE WGBS 流程 |
| OmicsWeave | 格式识别、转换为 bigBed（全部列）或 bigWig（一列） |

## 格式由来

bedMethyl 最早由 ENCODE 为 WGBS / RRBS 数据定义：在 BED9 后加两列，分别是覆盖该 CpG 的读段数和其中显示甲基化的百分比。ONT 的 `modkit` 沿用了这个名字和前 11 列的含义，又追加了 7 列计数，以区分修饰、未修饰、其他修饰、缺失、低置信度等调用结果，并在第 4 列写入修饰代码，使同一文件能同时包含 5mC、5hmC 等多种修饰。

两种方言前 11 列含义一致，主要差别在于 modkit 的附加计数和第 4、5 列的用法。

## 文件结构

### modkit pileup（18 列）

```text
chr1	10468	10469	m	12	+	10468	10469	255,0,0	12	75.00	9	3	0	0	1	0	2
chr1	10468	10469	h	12	+	10468	10469	255,0,0	12	8.33	1	3	8	0	1	0	2
```

| 列 | 名称 | 含义 |
|---|---|---|
| 1 | `chrom` | 参考序列名 |
| 2 | `start` | 位点起点，0-based |
| 3 | `end` | 位点右端，不包含 |
| 4 | `code` | 修饰代码，如 `m`（5mC）、`h`（5hmC）、`a`（6mA），或 ChEBI 编号；按 motif 输出时附带 motif 信息 |
| 5 | `score` | 等于 `Nvalid_cov` |
| 6 | `strand` | `+`、`-`；合并两条链（`--combine-strands`）时为 `.` |
| 7–8 | `thickStart` / `thickEnd` | 与 2、3 列相同 |
| 9 | `color` | 显示颜色，如 `255,0,0` |
| 10 | `Nvalid_cov` | 有效覆盖：`Nmod + Nother_mod + Ncanonical` |
| 11 | `percent_modified` | `Nmod / Nvalid_cov × 100`，范围 0–100 |
| 12 | `Nmod` | 调用为该修饰的读段数 |
| 13 | `Ncanonical` | 调用为未修饰碱基的读段数 |
| 14 | `Nother_mod` | 调用为同一碱基其他修饰的读段数 |
| 15 | `Ndelete` | 该位点为缺失的读段数 |
| 16 | `Nfail` | 调用概率低于阈值的读段数 |
| 17 | `Ndiff` | 该位点是其他碱基的读段数 |
| 18 | `Nnocall` | 是正确碱基但没有修饰调用的读段数 |

同一位点的不同修饰各占一行（上例 5mC 与 5hmC），因此同一坐标可能出现多行。较早版本的 modkit 用空格而不是 tab 分隔第 10–18 列，内容相同。

### ENCODE WGBS（11 列）

```text
chr1	10468	10469	.	14	+	10468	10469	255,0,0	14	79
```

| 列 | 名称 | 含义 |
|---|---|---|
| 1–3 | `chrom` / `start` / `end` | 位点坐标，0-based 半开 |
| 4 | `name` | 名称 |
| 5 | `score` | 覆盖度，上限截断为 1000 |
| 6 | `strand` | `+`、`-` 或 `.` |
| 7–8 | `thickStart` / `thickEnd` | 与 2、3 列相同 |
| 9 | `itemRgb` | 显示颜色，通常随甲基化水平变化 |
| 10 | `coverage` | 覆盖该位点的读段数 |
| 11 | `percent` | 显示甲基化的读段百分比，范围 0–100 |

## 关键概念

### 百分比不是 beta 值

两种方言的比例列都是 **0–100 的百分数**。许多下游模型使用 0–1 的 beta 值，需要自行除以 100。OmicsWeave 按原样保存数值，不做换算，也不从计数重新计算比例。

### 覆盖度决定可信度

同样是 50%，覆盖 2 条读段和覆盖 200 条读段的可信度完全不同。训练时常见的做法是同时提取比例和覆盖度两条轨道，再按覆盖度过滤或加权。这在 OmicsWeave 中对应两次 bedMethyl → bigWig 转换，各选一列。

### 链与 CpG

CpG 在两条链上是对称的：正链 C 与负链上相邻的 C（对应正链 G）通常甲基化状态一致。按链输出的文件中，同一个 CpG 会以相邻的两个位点出现在 `+` 和 `-` 上；`modkit --combine-strands` 则把它们合并为一个 `.` 链位点。

转为 bigWig 时，每个位置只能有一个值：如果两条链的记录在坐标上重叠，必须选择一条链。

### 多种修饰并存

modkit 文件可同时包含多种修饰代码，同一位点的 `m` 和 `h` 各占一行。转为 bigWig 时必须指定一种修饰；转为 bigBed 则全部保留。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 文件对象 | `BedMethyl`（纯文本）、`GzipBedMethyl`（gzip 或 bgzip） |
| 格式识别 | 每个数据行都是 BED9 加 modkit 的 9 列计数，或 ENCODE 的覆盖度与百分比；文件名不参与 |
| 方言 | modkit 18 列（含空格分隔计数的旧版写法）、ENCODE 11 列；按列数区分 |
| 区域查询 | 不直接支持；`query.open` 会报错并给出转换命令 |
| 转为 bigBed | 保留全部列，写入相应方言的 autoSql schema |
| 转为 bigWig | 选择一个数值列，可按修饰代码和链筛选 |
| 从 BAM 生成 | 不生成 bedMethyl；BAM 甲基化统计由 quant 输出胞嘧啶表和 bigWig |

## 保留全部字段：转为 bigBed

```sh
ow convert pileup.bed calls.bb
ow convert wgbs.bed.gz calls.bb --chrom-sizes genome.fa.fai
```

```python
from omicsweave.query import IntervalReader

with IntervalReader("calls.bb") as calls:
    print(calls.schema.field_count, calls.schema.defined_field_count)  # 18 9（modkit）
    print(calls.schema.auto_sql)                                       # 各列名称与说明
    hits = calls.fetch(("chr1", 10_000, 11_000))
    fields = [row.split("\t") for row in hits.annotation.to_list()]    # 第 4 列起
```

输出 bigBed 的前 9 列标记为 BED 标准字段，整份 autoSql 描述了每一列的名称、类型和含义，其他工具（如 UCSC 浏览器）也能据此显示各列。旧版 modkit 中空格分隔的计数在 bigBed 中统一以 tab 保存。

## 提取一条数值轨道：转为 bigWig

```sh
# 正链 5mC 比例（0–100）
ow convert pileup.bed 5mC.plus.bw --mod m --strand + --value percent

# 对应的有效覆盖度
ow convert pileup.bed 5mC.plus.coverage.bw --mod m --strand + --value coverage

# 合链输出的 CpG 文件
ow convert combined.bed 5mCpG.bw --mod m --strand .
```

```python
from omicsweave import convert

convert.run("pileup.bed", "5mC.plus.bw", mod="m", strand="+", value="percent")
convert.run("wgbs.bed.gz", "wgbs.coverage.bw", value="coverage")
```

| `value` | modkit 中的列 | ENCODE 中的列 |
|---|---|---|
| `percent`（默认，别名 `percent_modified`） | 11 `percent_modified` | 11 `percent` |
| `coverage`（别名 `valid_coverage`、`nvalid_cov`） | 10 `Nvalid_cov` | 10 `coverage` |
| `score` | 5，等于 `Nvalid_cov` | 5，截断到 1000 的覆盖度 |
| `nmod`、`ncanonical`、`nother_mod`、`ndelete`、`nfail`、`ndiff`、`nnocall` | 12–18 | 不适用，选择会报错 |

`value` 不区分大小写。筛选规则：

- `mod`：只保留第 4 列等于该代码的记录。文件中有多个修饰代码时必须指定；指定的代码不存在时报错，并列出文件中的代码。
- `strand`：只保留该链的记录，可选 `+`、`-`、`.`。`.` 只选择原本就标为合链的记录，**不会**替你把正负链合并。不指定时保留全部记录，但两条链的记录在坐标上重叠会报错，提示选择链。
- 筛选后仍有重叠（同一位置出现两个值）会报错，不会自动求和或求平均。

数值以 `float32` 原样保存。`chrom_sizes` 可选，接受两列长度表、`.fai` 或 `.2bit`；省略时按每条染色体最大的右端推断。

转换后读取：

```python
from omicsweave.query import TrackReader

with TrackReader("5mC.plus.bw") as meth:
    region = ("chr1", 10_000, 11_000)
    values = meth.signal(region).values          # 无调用的位置为 NaN
    beta = values / 100                          # 需要 0–1 的 beta 时自行换算
    bins = meth.bins(region, 100, ("mean", "coverage"))
```

`coverage` 统计量表示 bin 中有调用的位点比例，不是测序深度；测序深度来自 `--value coverage` 生成的另一条轨道。

## 输入要求

- 输入可以是纯文本、gzip 或 BGZF；`track`、`browser`、`#` 开头的行和空行跳过。
- 每行必须是 18 列或 11 列，两种方言**不能混在同一文件中**。
- 起点小于右端；第 4 列不能为空；链为 `+`、`-` 或 `.`；计数列为非负整数，比例列为数字。
- 不要求输入有序，转换按染色体（首次出现的顺序）和起点排序。

## 识别与显式指定

自动识别读取文件开头 64 KiB，要求其中每个数据行都符合某一方言的列结构，否则按普通 [BED](bed.md) 处理。识别结果不符合预期时，用格式对象显式指定；这样不符合方言的行会报出行号和具体原因，便于定位问题：

```python
from omicsweave.formats import BedMethyl

convert.run(BedMethyl("calls.txt"), "calls.bb")
```

CLI 对应 `--from bedmethyl`（压缩文件为 `--from gzip-bedmethyl`）。若希望把一个 bedMethyl 当作普通 BED 处理（不写 autoSql），可用 `--from bed`。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| `the items have N modification codes` | 有多种修饰；用 `--mod` 选择一种 |
| `overlaps ... on the other strand` | 两条链的记录重叠；用 `--strand +` 或 `--strand -` 分别输出 |
| `no item is on strand .` | 文件没有合链记录；`.` 不会自动合并正负链 |
| `value nmod is a modkit column` | ENCODE 方言没有 modkit 的计数列 |
| `11 fields; earlier records have 18` | 两种方言混在同一文件；分开转换 |
| 值在 0–100 之间 | 这是百分数；需要 beta 值时除以 100 |
| 输入还是 BAM | 先用 [quant 的甲基化功能](../quant.md) 从读段生成调用和统计 |

## 相关格式

- [BED](bed.md)：bedMethyl 的基础格式。
- [bigBed](bigbed.md)：保存 bedMethyl 全部列的可查询格式。
- [bigWig](bigwig.md)：保存其中一列数值信号的可查询格式。

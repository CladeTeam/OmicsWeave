# BED

[格式总览](README.md) · [查询](../query.md) · [转换](../convert.md) · [bigBed](bigbed.md) · [bedGraph](bedgraph.md)

BED 是描述**基因组区间**的制表文本：每行一个区间，前三列是染色体、起点和右端，后面可以跟名称、分值、链方向，乃至外显子结构。peaks、黑名单区域、增强子、CpG 岛、转录本模型，几乎所有“基因组上的一段段特征”都可以写成 BED，它也是生物信息学中流通最广的区间格式。

在 OmicsWeave 中，BED 是**转换源**：它没有索引，查一个区域要扫描整份文件，所以不直接用于区域查询，而是先用 `ow convert` 转为 [bigBed](bigbed.md)，再由 `IntervalReader` 按区域读取，或交给 `quant` 作为计数特征。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 基因组区间，及可选的名称、分值、链、blocks 等附加列 |
| 常见后缀 | `.bed`、`.bed.gz`；peaks 常见 `.narrowPeak`、`.broadPeak` |
| 存储 | 纯文本，可用 gzip / bgzip 压缩 |
| 坐标 | 0-based 半开 |
| 随机访问 | 文本本身不支持；转为 bigBed 后支持 |
| 外部索引 | 不需要（tabix 索引的 BED 也按普通文本读取） |
| 典型来源 | MACS2 / MACS3 peaks、ENCODE、UCSC Table Browser、bedtools、`Annotation.write_bed12` |
| OmicsWeave | 格式识别、转换为 bigBed、作为 consensus 的 mask 输入 |

## 格式由来

BED（Browser Extensible Data）最初是 UCSC 基因组浏览器用来上传自定义轨道的格式，因为简单，很快成为区间数据的通用交换格式。长期以来它只有 UCSC 网页上的非正式说明，各工具对细节（空格还是 tab、可否有 header、列数可否变化）理解不同；2022 年 GA4GH 发布了正式的 BED 规范，统一了大部分约定。

BED 的最大优点是“人人都能写”：`awk` 一行就能产出。代价是它没有索引、没有 schema，也没有染色体长度信息，这正是 bigBed 要补上的部分。

## 文件结构

### 列定义

前三列必需，后九列可选，但必须按顺序出现：要写第 6 列（strand），就必须先写第 4、5 列。

| 列 | 名称 | 含义 | 示例 |
|---|---|---|---|
| 1 | `chrom` | 染色体或参考序列名 | `chr1` |
| 2 | `chromStart` | 起点，0-based，包含 | `100` |
| 3 | `chromEnd` | 右端，不包含 | `200` |
| 4 | `name` | 特征名称 | `peak1` |
| 5 | `score` | 分值，UCSC 约定 0–1000 | `960` |
| 6 | `strand` | `+`、`-` 或 `.`（无链） | `+` |
| 7 | `thickStart` | 粗绘起点，基因模型中常为 CDS 起点 | `120` |
| 8 | `thickEnd` | 粗绘右端，常为 CDS 右端 | `480` |
| 9 | `itemRgb` | 显示颜色 `R,G,B`，或 `0` | `255,0,0` |
| 10 | `blockCount` | block（外显子）数 | `2` |
| 11 | `blockSizes` | 各 block 长度，逗号分隔 | `100,150,` |
| 12 | `blockStarts` | 各 block 相对 `chromStart` 的起点 | `0,250,` |

按列数习惯称为 BED3、BED6、BED12 等：

```text
# BED3：只有区间
chr1	100	200

# BED6：加上名称、分值、链
chr1	100	200	peak1	960	+
chr1	300	450	peak2	512	-

# BED12：一个两外显子的转录本
chr1	100	500	tx1	0	+	120	480	0	2	100,150,	0,250,
```

### BED12 的 blocks

BED12 用 blocks 在一行内描述“拼接”的特征，最常见的就是转录本：

```text
chromStart=100                                         chromEnd=500
   │                                                        │
   ▼                                                        ▼
   [=== block 1: 100–200 ===]------ 内含子 ------[=== block 2: 350–500 ===]
   blockStarts[0]=0, blockSizes[0]=100           blockStarts[1]=250, blockSizes[1]=150
```

每个 block 的绝对坐标是 `chromStart + blockStarts[i]` 到 `chromStart + blockStarts[i] + blockSizes[i]`。第一个 block 必须从 0 开始，最后一个 block 必须到 `chromEnd` 结束。

### BEDn+m：扩展列

在标准列之后再加自定义列的写法称为 BEDn+m：前 n 列遵循 BED 标准，后 m 列由各自格式定义。常见例子：

| 格式 | 结构 | 附加列 |
|---|---|---|
| ENCODE narrowPeak | BED6+4 | `signalValue`、`pValue`、`qValue`、`peak`（峰顶相对起点的偏移） |
| ENCODE broadPeak | BED6+3 | `signalValue`、`pValue`、`qValue` |
| [bedMethyl](bedmethyl.md) | BED9+2 / BED9+9 | 覆盖度、甲基化比例或修饰计数 |

### header 行

数据行之前可以有：

- `track ...`：UCSC 浏览器的轨道设置，如 `track name=peaks color=255,0,0`。
- `browser ...`：浏览器显示设置。
- `#` 开头的注释。

## 关键概念

### 0-based 半开坐标

BED 是 0-based 半开坐标：`chr1 100 200` 表示 `[100, 200)`，覆盖第 100–199 号碱基（从 0 数），长度正好是 `end - start = 100`。

容易混淆的地方在于，基因组浏览器界面、GTF / GFF3、VCF、SAM 都使用 1-based 闭区间。同一段区域：

| 表示方式 | 写法 |
|---|---|
| BED、OmicsWeave 查询 key | `chr1 100 200` |
| UCSC 浏览器地址栏、GTF | `chr1:101-200` |

从 1-based 闭区间换成 BED 时，起点减一、右端不变。OmicsWeave 的所有查询都使用 0-based 半开坐标，与 BED 一致，因此 BED 中的区间可以直接作为 `Region` 使用。

### 区间可以重叠

BED 不要求区间互不重叠，也不要求排序。peaks、转录本、重复元件经常彼此重叠。这一点和 [bedGraph](bedgraph.md) 不同：bedGraph 表示信号，每个位置只能有一个值。

### 没有染色体长度

BED 文件本身不记录染色体有多长。转换为 bigBed 时需要染色体长度：可以由 `chrom_sizes` 提供，也可以退而使用每条染色体观测到的最大右端。

### BED 与基因注释

BED12 能表示转录本的外显子结构，但不记录基因与转录本的从属关系、基因名、biotype 等。需要这些信息时使用 [GTF](gtf.md) / [GFF3](gff3.md) 的 `Annotation`；反过来，`Annotation.write_bed12` 可以把转录本导出为 BED12。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 文件对象 | `Bed`（纯文本）、`GzipBed`（gzip 或 bgzip） |
| 格式识别 | `formats.detect` 按内容判断；文件名不参与 |
| 区域查询 | 不直接支持；`query.open` 和 `IntervalReader` 会报错并给出转换命令 |
| 转换 | BED → bigBed，保留全部附加列 |
| 导出 | `Annotation.write_bed12` 把转录本写成 BED12 |
| 其他用途 | `ConsensusOptions(mask=...)` 接受 BED 文件作为屏蔽区域 |
| 定量 | `quant.Peaks` 需要 bigBed，BED 先转换 |

## 转为 bigBed

```sh
ow convert peaks.bed peaks.bb --chrom-sizes genome.fa.fai
ow convert peaks.bed.gz peaks.bb --chrom-sizes genome.2bit
ow convert peaks.narrowPeak peaks.bb --from bed --chrom-sizes genome.chrom.sizes
```

```python
from omicsweave import convert

report = convert.run("peaks.bed", "peaks.bb", chrom_sizes="genome.chrom.sizes")
print(report.target, report.records)
```

转换规则：

- **输入**：纯文本、gzip、BGZF 均可。`track`、`browser`、`#` 开头的行和空行跳过。
- **分隔符**：含 tab 的行按 tab 切分；不含 tab 的行按连续空白切分。名称中含空格时应使用 tab 分隔。
- **列数**：至少三列，且所有记录的列数必须相同；第 4 列起原样保存，以 tab 连接。
- **坐标**：起点、右端须是 `uint32` 范围内的非负整数，`start > end` 报错。
- **排序**：不要求输入有序。转换按染色体（首次出现的顺序）和起点排序，起点相同时保持输入顺序；数据超过内存预算时使用临时文件外排序。
- **重叠**：允许。
- **染色体长度**：`chrom_sizes` 接受两列长度表、`.fai` 或 `.2bit`。提供后，记录中出现但长度表中没有的染色体会报错，右端超过染色体长度也会报错，不会自动裁剪。省略时以每条染色体最大的右端作为长度。输出只包含实际有记录的染色体。
- **schema**：前 `min(列数, 12)` 列标记为 BED 标准字段，其余为扩展字段；BED 转换不写 autoSql。`IntervalReader.schema.fields` 的标准回退名称把第 9、12 列分别称为 `reserved`、`chromStarts`（对应上表的 `itemRgb`、`blockStarts`）。
- 可选 `--sidecar` 同时写出 [sidecar](sidecar.md)，普通查询不需要。

转换后查询：

```python
from omicsweave.query import IntervalReader

with IntervalReader("peaks.bb") as peaks:
    hits = peaks.fetch(("chr1", 150, 350))
    print(hits.start, hits.end)          # 与区域重叠的完整记录，不裁剪
    print(hits.annotation.to_list())     # ['peak1\t960\t+', 'peak2\t512\t-']
```

`annotation` 中是第 4 列起的原始文本，不会自动拆分成 name、score、strand。需要某一列时自行 `split("\t")`。`fetch(..., clip=True)` 可把返回坐标截到查询区域，附加字段保持原样。

## 识别与显式指定

自动识别只检查开头最多 64 KiB 中的完整行：样本中有数据行，且每个数据行都恰好四列、起点和右端可解析为 `uint32`、最后一列可解析为浮点数时，会识别为 [bedGraph](bedgraph.md)；`track type=bedGraph` 也可直接指定该类型。例如 `chr1 100 200 5` 既可能是带分值的区间，也可能是信号。区间语义更合适时，显式指定格式：

```python
from omicsweave import convert
from omicsweave.formats import Bed, detect

print(detect("regions.txt"))   # 可能是 BedGraph('regions.txt')
convert.run(Bed("regions.txt"), "regions.bb", chrom_sizes="genome.chrom.sizes")
```

CLI 对应 `--from bed`。符合 [bedMethyl](bedmethyl.md) 列结构的文件会被识别为 bedMethyl，转换时额外写入相应 autoSql；同样可以用 `--from bed` 按普通 BED 处理。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| `bed is not read by region` | BED 不直接查询；按提示用 `ow convert` 转为 bigBed |
| 被识别为 bedGraph | 探测样本符合四列 bedGraph 条件，或有相应 track 声明；用 `Bed(path)` 或 `--from bed` 显式指定 |
| `N fields; earlier records have M` | 所有行列数必须一致；先补齐或裁掉多余列 |
| `chromosome ... is not in the chromosome sizes` | 长度表与 BED 的染色体命名不一致（如 `chr1` 与 `1`），或长度表不完整 |
| `end ... past the end of ...` | 区间超出染色体长度；检查参考版本，或先自行裁剪 |
| 查询结果坐标超出查询区域 | `IntervalReader` 默认返回完整记录；需要裁剪时使用 `fetch(region, clip=True)` |
| 与浏览器显示差 1 | 浏览器显示 1-based；BED 起点比浏览器显示的小 1 |

## 相关格式

- [bigBed](bigbed.md)：BED 的带索引二进制形式，OmicsWeave 查询区间的格式。
- [bedGraph](bedgraph.md)：形似 BED4 的数值信号格式。
- [bedMethyl](bedmethyl.md)：BED9 加甲基化计数的扩展格式。
- [GTF](gtf.md)、[GFF3](gff3.md)：带基因层级的注释格式，可导出 BED12。

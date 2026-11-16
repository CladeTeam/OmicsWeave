# bigWig

[格式总览](README.md) · [查询](../query.md) · [转换](../convert.md) · [定量](../quant.md) · [bedGraph](bedgraph.md) · [wig](wig.md)

bigWig 是保存**全基因组连续数值信号**的二进制格式：覆盖度、ChIP / ATAC 富集、RNA-seq 链特异信号、保守性分数、甲基化比例等都常以 bigWig 发布。它把信号压缩成块，并自带索引和多层预汇总，因此可以只读取需要的一小段，也能快速得到宽窗口的统计值。

在 OmicsWeave 中，bigWig 是数值信号的标准查询格式：bedGraph、wig 转换为 bigWig 后再按区域读取。bedMethyl 的区域查询转换提示指向 bigBed，以保留全部列；要提取一列数值，可显式选择 bigWig 目标。`quant` 计算的覆盖度和甲基化轨道也直接输出 bigWig。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 基因组坐标上的浮点信号，每个区间一个值 |
| 常见后缀 | `.bw`、`.bigwig`、`.bigWig` |
| 存储 | 二进制 BBI 容器，数据块用 zlib 压缩 |
| 坐标 | 0-based 半开，与 bedGraph 相同 |
| 随机访问 | 支持，文件内置 R 树索引 |
| 外部索引 | 不需要 |
| 典型来源 | deepTools `bamCoverage`、UCSC `bedGraphToBigWig` / `wigToBigWig`、ENCODE、`ow quant --coverage` |
| OmicsWeave | 区域查询、分箱统计、直接写入、作为转换目标和 quant 输出 |

## 格式由来

bigWig 由 UCSC 基因组浏览器团队设计（Kent et al., *Bioinformatics* 2010），与 [bigBed](bigbed.md) 共用同一种称为 **BBI**（Big Binary Indexed）的容器。设计目标是让浏览器通过 HTTP Range 请求只取屏幕上可见的部分，所以它天然适合“给定区域，读一段信号”的访问模式，这也正是训练数据加载需要的。

它的文本前身是 [bedGraph](bedgraph.md) 和 [wig](wig.md)。文本格式便于生成和阅读，但查询一个区域需要从头扫描；bigWig 保存相同的“区间 → 数值”信息，只是额外加入了压缩、索引和汇总层。

## 文件结构

OmicsWeave writer 生成的 bigWig 从前到后包含以下段；读取其他文件时以 header 的偏移为准：

```text
┌──────────────────────────────┐
│ 固定 header（magic、版本、各段偏移）│
│ zoom header 表及预留空间     │  每个汇总层的分辨率与位置
│ 全文件 summary               │  覆盖碱基数、最小 / 最大值、总和、平方和
│ extension header             │
│ 染色体 B+ 树                 │  染色体名 → 编号、长度
│ 数据 section 数（uint64）    │  OmicsWeave 每个原始数据块一个 section
│ 数据块（zlib 压缩）          │  原始信号
│ R 树索引                     │  区域 → 数据块位置
│ zoom 层 1 数据 + 索引        │
│ zoom 层 2 数据 + 索引        │  分辨率逐层变粗
│ ...                          │
│ 尾部 magic                   │
└──────────────────────────────┘
```

**数据块**中的条目有三种写法，对应文本格式的三种段：

| 块类型 | 每个条目 | 来源 |
|---|---|---|
| bedGraph | `start, end, value` | bedGraph 行、wig 中的 bedGraph 行 |
| variableStep | `start, value`，跨度固定 | wig `variableStep` |
| fixedStep | `value`，起点和步长固定 | wig `fixedStep` |

三者描述的都是“某段区间上的一个数值”，读取时没有区别。

**zoom 层**是对原始数据的预汇总：每条记录覆盖最多 `reduction` 个碱基（例如 830 bp、3320 bp……），染色体末端的记录可以更短，保存这段内的有效碱基数、最小值、最大值、总和和平方和。查询一个 1 Mb 宽的 bin 时，可以直接读几十条 zoom 记录，而不必解压几百个原始数据块。默认推导的层只在记录数逐层减少时保留，空文件和部分小文件没有 zoom 层。

**R 树索引**把基因组区间映射到数据块的文件位置。查询时先在索引中找出与区域重叠的块，再只读取、解压这些块。

## 关键概念

### 没有数据 ≠ 0

bigWig 只存有值的区间，未列出的位置就是“没有数据”。这和“信号为 0”不同：覆盖度轨道里没有读段的地方通常不写，但保守性分数缺失的地方可能确实没有测量。OmicsWeave 用 NaN 表示没有数据，是否把它视为 0 由使用者按数据含义决定：

```python
import numpy as np

values = track.signal(region).values
values = np.nan_to_num(values, nan=0.0)  # 仅当“缺失即 0”符合数据含义时
```

### 原始值与统计值

同一个文件可以回答两类问题：

- **逐碱基值**：每个位置一个数，用于碱基分辨率模型。
- **分箱统计**：把区域切成固定宽度的 bin，每个 bin 给出 mean、max、coverage 等，用于较粗分辨率的模型输入或可视化。

分箱时，mean 只对**有数据的碱基**求均值。一个 100 bp 的 bin 中只有 20 bp 的值为 5，则 `mean=5`、`sum=100`、`coverage=0.2`；如果希望把缺失视为 0，应使用 `sum / bin 宽度 = 1`。

### 精确统计与 zoom 近似

精确统计读取原始数据块；近似统计选择 reduction 不超过半个 bin 宽度的最粗 zoom 层，并按 zoom 记录落在 bin 内的比例加权。宽 bin 下近似统计读取的数据少得多，但结果是估计值。

### 链方向

bigWig 本身没有链的概念。链特异数据（RNA-seq、CAGE、PRO-seq）通常拆成正链、负链两个文件，有的工具会把负链信号写成负数。OmicsWeave 用 `StrandedTrackReader` 按 key 的链选择文件，数值按文件原样返回。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 格式识别 | `formats.detect` 按 magic 识别为 `BigWig` |
| 逐碱基信号 | `TrackReader.signal`，`float32`，无数据为 NaN |
| 分箱统计 | `TrackReader.bins`：mean、min、max、sum、std、coverage，可一次请求多个 |
| 原始区间 | `TrackReader.intervals`：文件中存储的 runs，默认裁剪到查询区域，`clip=False` 保留原始起止 |
| 多段 key | `signal` / `bins` 支持 `Spliced`，适合转录本信号 |
| 链方向 | 负链 key 按 5′ → 3′ 反转顺序；分链文件用 `StrandedTrackReader` |
| 文件信息 | `chroms`、`summary`、`zoom_reductions` |
| 缓存与并发 | 可接 `BlockCache`；解码释放 GIL；可 `clone()`、pickle |
| 写入 | `BigWigWriter`，可多线程压缩块，可选同时写 sidecar |
| 转换生成 | bedGraph、wig、bedMethyl（选一个数值列）→ bigWig |
| quant 输出 | 覆盖度、甲基化 beta / coverage 轨道 |

限制：

- 只读本地文件，不直接读取 HTTP URL。
- 只支持 little-endian 文件（常见工具生成的文件都是 little-endian）。
- 写出 BBI version 4，与 UCSC 工具兼容。
- 当前不提供 bigWig → bedGraph / wig 的文本导出。

## 查询

### 逐碱基信号

```python
from omicsweave.query import TrackReader

with TrackReader("signal.bw") as track:
    print(dict(track.chroms))                 # {'chr1': 248956422, ...}
    signal = track.signal(("chr1", 0, 1000))  # float32 [1000]
    print(signal.values[:10])
```

区域可以超过染色体末端，超出部分为 NaN；染色体名必须与文件中的一致（`chr1` 和 `1` 不会自动互换）。

### 分箱统计

```python
with TrackReader("signal.bw") as track:
    bins = track.bins(("chr1", 0, 1000), 100, ("mean", "max", "coverage"))
    print(bins.values.shape)  # (3, 10)：[统计量数, bin 数]
    print(bins["mean"])       # 按统计量名取一行
    print(bins.edges)         # [0, 100, ..., 1000]，相对 key 起点的局部坐标
```

| 统计量 | 含义 | 空 bin |
|---|---|---|
| `mean` | 有数据碱基的均值 | NaN |
| `min` / `max` | 有数据碱基的最小值 / 最大值 | NaN |
| `sum` | 按覆盖碱基数加权的总和 | NaN |
| `std` | 有数据碱基的样本标准差；原始数据路径只有一个有效碱基时为 0，zoom 路径估计的有效碱基数不大于 1 时为 NaN | NaN |
| `coverage` | 有数据碱基占 bin 的比例 | 0 |

常用选项：

- `phase`：分箱起点偏移。默认从 key 的局部位置 0 开始；对齐全基因组网格时可设 `phase=region.start % resolution`。
- `exact=False`：使用 zoom 层近似，适合宽 bin；没有合适的 zoom 层时自动读原始数据。`exact=True` 强制使用原始数据，因此单碱基 `std=0` 的规则也适用于回退到原始数据的查询。
- `out=`：传入 `float32 [统计量数, bin 数]` 的 C-contiguous 数组，结果直接写入，适合拼 batch。

当前实现的精确统计从 full-resolution 数据累积；近似统计按上述 reduction 规则选择 zoom 层。结果精度还取决于文件的 zoom 布局和 bin 边界，不应把近似统计当作精确值。

### 原始区间

```python
with TrackReader("signal.bw") as track:
    runs = track.intervals(("chr1", 120, 250))
    print(runs.start, runs.end, runs.value)  # uint32、uint32、float32
    stored = track.intervals(("chr1", 120, 250), clip=False)
```

返回文件中存储的“区间 → 值”记录，并始终保持基因组顺序。`clip=True` 是默认值，起止坐标裁剪到查询区域；`clip=False` 保留 run 在文件中的原始起止，`value` 不变。两种模式都只返回重叠的 runs。适合稀疏信号，或需要自行做加权汇总时使用。

`intervals` 只接受连续 `Region`，负链不会反转区间列表；`signal` / `bins` 则支持 `Spliced` 并按链定向。`StrandedTrackReader.intervals` 同样接受 `clip`，并按 key 的链选择文件。

### 转录本与链方向

```python
from omicsweave.query import StrandedTrackReader
from omicsweave.types import Spliced

transcript = Spliced("chr1", [(1000, 1200), (1500, 1700)], strand="-")

with StrandedTrackReader("rna.plus.bw", "rna.minus.bw") as rna:
    values = rna.signal(transcript).values  # 400 个值，按负链 5′ → 3′
```

`StrandedTrackReader` 要求 key 带有明确的 `+` 或 `-`；`.` 无法决定读取哪个文件，会报错。单文件 `TrackReader` 也接受负链 key，它只反转输出顺序。

### 文件信息

```python
with TrackReader("signal.bw") as track:
    print(track.summary)          # BigWigSummary(...)；文件未存总汇总时为 None
    print(track.zoom_reductions)  # 各层 reduction（每条记录的最大跨度），由细到粗
```

`summary` 由文件 header 指向的总汇总读取；未提供时为 `None`。存在且 `bases_covered > 0` 时，可用 `sum / bases_covered` 得到全基因组有数据碱基的均值，例如用于训练前的信号标准化。

## 生成 bigWig

### 从文本信号转换

```sh
ow convert signal.bedGraph signal.bw --chrom-sizes genome.fa.fai
ow convert signal.wig signal.bw --chrom-sizes genome.2bit
ow convert pileup.bed 5mC.plus.bw --mod m --strand + --value percent
```

`--chrom-sizes` 接受两列长度表、`.fai` 或 `.2bit`；省略时以每条染色体观测到的最大右端作为长度。转换会先排序输入，但区间不能重叠，也不能超出声明的染色体长度。各来源的细节见 [bedGraph](bedgraph.md)、[wig](wig.md)、[bedMethyl](bedmethyl.md)。

### 从 BAM 计算覆盖度

```sh
ow quant reads.bam -o results --coverage
ow quant reads.bam -o results --coverage --strands split --library 2 --scale cpm
```

BAM → bigWig 是计算覆盖度，属于 [quant](../quant.md)，不是格式转换：需要决定计数单位、过滤条件、是否分链和缩放方式。

### 从内存写入

```python
import numpy as np
from omicsweave.formats import BbiOptions, BigWigWriter

chroms = {"chr1": 1000, "chr2": 500}
with BigWigWriter("example.bw", chroms, sidecar=True, threads=4) as writer:
    writer.add("chr1", [100, 200], [150, 300], [2.5, 8.0])   # 数组形式
    writer.add("chr2", 0, 100, 1.0)                          # 标量形式
```

- `chroms` 的顺序就是文件中的染色体顺序；记录必须按这个顺序、按起点排序，且互不重叠。
- writer 用上一条记录的右端检查下一条起点，同一染色体上 `start < 上一条 end` 会报排序错误；不会合并重叠记录。
- `add` 的 start、end、value 可以是标量或等长数组，会自动广播；整段写入 NumPy 数组比逐条调用快得多。
- `options=BbiOptions(...)` 可调整块大小（`items_per_slot`，默认 1024）、索引扇出（`block_size`，默认 256）、压缩等级（`level`，0–12，默认 6）、zoom 层数（`zoom_levels`，默认 10）或显式指定 `zoom_reductions`。日常使用默认值即可。
- `threads` 大于 1 时在工作线程上压缩数据块，输出文件与线程数无关。
- 离开 `with` 时完成文件；块内发生异常时放弃这个未完成的文件。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 信号中有大量 NaN | 这些位置没有数据；按数据含义决定是否填 0 |
| mean 比预期高 | mean 只对有数据碱基求均值；缺失视为 0 时用 `sum / 宽度` |
| 找不到染色体 | 检查 `track.chroms`，使用文件中的实际名称 |
| 宽窗口统计较慢 | 使用 `exact=False` 读取 zoom 层 |
| 负链信号是负数 | 文件本身就这样存储；OmicsWeave 不改变数值符号 |
| 转换时报区间重叠 | bigWig 每个位置只能有一个值；先合并或拆分源数据 |

## 相关格式

- [bedGraph](bedgraph.md)、[wig](wig.md)：bigWig 的文本前身，转换为 bigWig 后查询。
- [bigBed](bigbed.md)：同一 BBI 容器中的区间记录，适合 peaks 和注释特征。
- [bedMethyl](bedmethyl.md)：可提取一个数值列成为 bigWig。
- [Sidecar](sidecar.md)：bigWig 数据块布局的可选伴随文件，普通查询不需要。

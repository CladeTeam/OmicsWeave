# .hic

[格式总览](README.md) · [查询](../query.md) · [转换](../convert.md) · [定量](../quant.md) · [cool](cool.md) · [mcool](mcool.md) · [4DN pairs](pairs.md) · [Juicer short](juicer-short.md)

`.hic` 是 Aiden 实验室 Juicer 工具链用来保存 **Hi-C / Micro-C 染色质接触矩阵**的二进制格式。一个文件包含全部染色体对、多种分辨率的接触矩阵，按块压缩并带有索引，因此可以只取出某个二维窗口，例如“chr1 上 0–1 Mb 对 0–1 Mb，5 kb 分辨率”的一块矩阵。

在 OmicsWeave 中，`.hic` 是接触数据的标准查询格式：[cool](cool.md)、[mcool](mcool.md)、[4DN pairs](pairs.md) 和 [Juicer short](juicer-short.md) 都转换为 `.hic` 后按窗口读取，`quant` 从 Hi-C BAM 计算的接触也直接写成 `.hic`。查询返回原始计数，可以直接作为模型的二维输入。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 染色体对之间分箱后的接触计数 |
| 常见后缀 | `.hic` |
| 存储 | 二进制，矩阵切成块，每块 zlib 压缩 |
| 坐标 | BP 分辨率下为碱基坐标；FRAG 分辨率下为酶切片段编号 |
| 随机访问 | 支持，文件 footer 中有主索引和各矩阵的块索引 |
| 外部索引 | 不需要 |
| 版本 | OmicsWeave 读取 v8 / v9，写出 v9（默认） |
| 典型来源 | Juicer `pre`、4DN / ENCODE 数据门户、`ow convert`、`ow quant --hic` |
| OmicsWeave | 窗口查询（稠密矩阵、非零像素）、直接写入、作为转换目标和 quant 输出 |

## 格式由来

Hi-C 实验测量基因组中任意两个位点在三维空间里相互靠近的频率。比对后的原始结果是一组“接触对”（每对读段两端的位置），常以 [4DN pairs](pairs.md) 或 [Juicer short](juicer-short.md) 文本保存。把基因组切成等宽 bin，统计每对 bin 之间的接触数，就得到一个对称的接触矩阵。

全基因组矩阵非常大：人类基因组 5 kb 分辨率约有 60 万个 bin，矩阵有 3.6 × 10¹¹ 个格子。Juicer 团队设计 `.hic`（Durand et al., *Cell Systems* 2016）来保存这种稀疏矩阵：只存非零格子，切块压缩，并预先计算从 2.5 Mb 到 5 kb 的多种分辨率，供 Juicebox 浏览器平滑缩放。另一种常用存储是基于 HDF5 的 [cool / mcool](cool.md)。

## 文件结构

```text
┌─────────────────────────────────────┐
│ header                              │  magic "HIC"、版本、footer 位置、genome、
│                                     │  属性、染色体表、BP / FRAG 分辨率、酶切位点
├─────────────────────────────────────┤
│ body：每个染色体对的矩阵            │
│   matrix (chr1, chr1)               │
│     resolution BP 2.5 Mb → blocks   │
│     resolution BP 1 Mb   → blocks   │
│     ...                             │
│     resolution BP 5 kb   → blocks   │  每块 zlib 压缩，记录 (x bin, y bin, 计数)
│   matrix (chr1, chr2)               │
│   ...                               │
├─────────────────────────────────────┤
│ footer                              │  主索引：染色体对 → 矩阵位置
│                                     │  expected 值、归一化向量（可选）
└─────────────────────────────────────┘
```

层级关系是 **矩阵 → 分辨率 → 块**：

- **矩阵**：每个染色体对一个，例如 `(chr1, chr1)` 是染色体内（cis / intra），`(chr1, chr2)` 是染色体间（trans / inter）。
- **分辨率**：同一矩阵按多种 bin 宽度保存。单位可以是 **BP**（碱基对）或 **FRAG**（限制性酶切片段，bin 编号按片段计数）。
- **块**：每个分辨率的矩阵被切成方块网格，每块单独压缩。块索引记录每块在文件中的位置和大小，查询时只读与窗口重叠的块。

### 对称矩阵只存一半

染色体内的接触矩阵是对称的：`(i, j)` 和 `(j, i)` 是同一组接触。`.hic` 只保存其中一半，读取时由读取方补全。OmicsWeave 的查询总是返回完整的窗口，无论文件存的是哪一半。

### v8 与 v9

| | v8 | v9 |
|---|---|---|
| 染色体长度 | 32-bit | 64-bit |
| footer 长度字段 | 32-bit | 64-bit |
| 归一化向量索引 | 在 footer 中 | header 额外记录其位置和长度 |
| 染色体内块布局 | 方形网格 | 可按到对角线的距离分层排列（`v9-depth-base` 属性） |

v9 的对角线布局让沿对角线滑动的窗口读取更少的块，这正是 Hi-C 模型常见的访问方式。v8 之前的版本以及更新的实验版本，OmicsWeave 不读取。

### 归一化与 expected

Hi-C 原始计数存在测序深度、GC、可比对性等偏差，常用 KR、VC、VC_SQRT、SCALE 等方法计算每个 bin 的校正向量；expected 是按距离对角线的远近统计的平均接触值，用于计算 observed / expected。`.hic` 可以在 footer 中保存这些向量，但它们是**分析结果**，不是矩阵本身。OmicsWeave 只读取原始计数（Juicer 中 `NONE` 归一化），不读取也不计算归一化向量和 expected；写出的文件中这两部分为空列表。

## 关键概念

### 分辨率与 bin

BP 分辨率 `r` 下，位置 `p` 属于第 `p // r` 个 bin。一个窗口包含与区域重叠的**所有** bin，所以区域边界不需要与 bin 对齐：

```python
from omicsweave.types import ContactWindow

window = ContactWindow(("chr1", 2_500, 1_002_500), ("chr1", 0, 1_000_000), 5_000)
print(window.shape)     # (201, 200)：行区域跨越 201 个 bin
print(window.row_bins)  # range(0, 201)
```

分辨率必须是文件中已有的一种，查询不会临时重新分箱。

### 稠密矩阵与稀疏像素

同一个窗口有两种取法：

- **稠密矩阵**：`[行 bin 数, 列 bin 数]` 的 `float32` 数组，没有接触的格子为 0。适合固定大小的模型输入。
- **非零像素**：`row`、`col`、`value` 三列，`row` / `col` 是 bin 起点坐标（窗口单位下的碱基或片段编号），不是矩阵下标。适合稀疏数据或大窗口。

### 方向

稠密矩阵的行、列各自按对应区域的链方向排列：行区域为 `-` 时行轴反转，列区域同理。这样可以按转录方向对齐接触图。像素始终保持基因组坐标。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 格式识别 | `formats.detect` 按 magic 识别为 `Hic` |
| 稠密矩阵 | `ContactReader.matrix`，`float32 [行, 列]`，可传 `out=` |
| 非零像素 | `ContactReader.pixels`，按 (row, col) 排序 |
| 单位 | BP；文件含 FRAG 分辨率时也可按片段查询 |
| 染色体间 | 行、列可来自不同染色体 |
| 文件信息 | `version`、`genome`、`attributes`、`chroms`、`maps`、`resolutions()`、`sites()`、`fragment_count()` |
| 缓存与并发 | 可接 `BlockCache`；解码释放 GIL；可 `clone()`、pickle |
| 写入 | `HicBuilder`（从接触位置分箱）、`HicWriter`（从已分箱像素），可多线程压缩、可选 sidecar |
| 转换生成 | cool、mcool、4DN pairs、Juicer short → `.hic` |
| quant 输出 | `ow quant --hic`：从 Hi-C / Micro-C BAM 生成 pairs 和 `.hic` |

限制：

- 只返回原始计数，不应用 KR / VC / SCALE 等归一化，不提供 expected。
- convert 和 quant 生成的文件只有 BP 分辨率，不生成 FRAG 分辨率，也不包含 Juicer 的 `All` 伪染色体。
- 当前不提供 `.hic` → cool / mcool 的反向导出。

## 查询

### 打开文件并查看内容

```python
from omicsweave.query import ContactReader

with ContactReader("contacts.hic") as contacts:
    print(contacts.version, contacts.genome)  # 9 hg38
    print(contacts.chroms)                    # {'chr1': 248956422, ...}
    print(contacts.resolutions())             # (2500000, 1000000, ..., 5000)
    print(contacts.maps[:3])                  # (('chr1', 'chr1'), ('chr1', 'chr2'), ...)
    print(contacts.attributes)                # {'software': ...}
```

`resolutions()` 默认列出 BP 分辨率，顺序与文件一致（Juicer 工具按由粗到细）。`maps` 是有接触矩阵的染色体对；对其他染色体对查询会得到全 0。

### 稠密矩阵

```python
from omicsweave.types import ContactWindow

with ContactReader("contacts.hic") as contacts:
    window = ContactWindow(("chr1", 0, 1_000_000), ("chr1", 0, 1_000_000), 5_000)
    matrix = contacts.matrix(window)
    print(matrix.values.shape)  # (200, 200)
    print(matrix.values.dtype)  # float32
```

窗口在对角线上时返回完整的对称矩阵。染色体间窗口：

```python
trans = ContactWindow(("chr1", 0, 2_000_000), ("chr2", 0, 1_000_000), 25_000)
block = contacts.matrix(trans)  # (80, 40)
```

拼 batch 时可传入预分配的 `float32`、C-contiguous、形状等于 `window.shape` 的数组：

```python
import numpy as np

buffer = np.empty(window.shape, dtype=np.float32)
contacts.matrix(window, out=buffer)
```

### 非零像素

```python
with ContactReader("contacts.hic") as contacts:
    pixels = contacts.pixels(window)
    print(len(pixels))
    print(pixels.row, pixels.col, pixels.value)  # uint32, uint32, float32
```

染色体内窗口中的对称格子会以两个像素出现（`(i, j)` 和 `(j, i)`）。转成窗口内下标：

```python
row_index = (pixels.row - window.row_bins.start * window.resolution) // window.resolution
```

### 按链方向读取

```python
window = ContactWindow(("chr1", 0, 1_000_000, "-"), ("chr1", 0, 1_000_000), 5_000)
matrix = contacts.matrix(window)  # 行轴反转，列轴不变
```

### FRAG 分辨率

已有文件含酶切片段分辨率时：

```python
from omicsweave.types import BinUnit

with ContactReader("contacts.hic") as contacts:
    print(contacts.resolutions("FRAG"))   # 例如 (500, 200, 100, 50, 20, 5, 2, 1)
    print(contacts.fragment_count("chr1"))
    sites = contacts.sites("chr1")         # uint32 酶切位点位置
    window = ContactWindow(("chr1", 0, 1000), ("chr1", 0, 1000), 1, unit=BinUnit.FRAG)
    matrix = contacts.matrix(window)
```

FRAG 窗口中的区域坐标是片段编号，不是碱基位置。

## 生成 .hic

### 从其他接触格式转换

```sh
ow convert sample.cool contacts.hic -t 4 -m 2G
ow convert sample.mcool contacts.hic --resolutions 5000,25000
ow convert sample.pairs.gz contacts.hic --resolutions 5000,25000
ow convert sample.short contacts.hic --from juicer-short \
  --chrom-sizes genome.chrom.sizes --resolutions 5000,25000
```

- cool / mcool 保留源文件已有的分辨率和原始 count，不读取 balance weight。
- pairs / Juicer short 按指定分辨率分箱；默认是 Juicer `pre` 的九档 2.5 Mb、1 Mb、500 kb、250 kb、100 kb、50 kb、25 kb、10 kb、5 kb。较粗分辨率必须是最细分辨率的整数倍。

细节见各源格式文档。

### 从 BAM 计算

```sh
ow quant hic.bam -o hic-results --hic --hic-resolutions 5000,25000 --genome hg38
```

quant 从两端 mate 的 5′ 位置得到接触，做 MAPQ 过滤和精确去重，输出 `hic/contacts.pairs.gz` 和 `hic/contacts.hic`。过滤规则见 [quant](../quant.md)。

### 从内存写入

已有接触位置（每对一条）时用 `HicBuilder`，关闭时按 header 中的每个分辨率分箱，与 Juicer `pre` 的做法一致：

```python
import numpy as np
from omicsweave.formats import HicBuilder, HicHeader

header = HicHeader("hg38", {"chr1": 1_000_000, "chr2": 500_000}, bp_resolutions=(25_000, 5_000))
with HicBuilder("example.hic", header, threads=4) as writer:
    writer.add("chr1", 10_000, "chr1", 50_000)                # 一条接触，count 默认 1
    writer.add("chr1", np.array([100, 200]), "chr2", np.array([300, 400]), count=2.0)
```

`add` 的所有参数可以是标量或数组并自动广播；染色体可以是名称，也可以是 header 中的编号数组。

已有分箱像素时用 `HicWriter.write_matrix`，x / y 是 bin 编号：

```python
from omicsweave.formats import HicWriter

with HicWriter("binned.hic", header) as writer:
    writer.write_matrix("chr1", "chr1", 5_000, x=[0, 1, 2], y=[0, 3, 2], value=[5.0, 1.0, 2.0])
```

格子可以任意顺序、任意方向（上三角或下三角）给出，落在同一存储位置的值会相加。

`HicHeader` 的字段：

| 字段 | 含义 |
|---|---|
| `genome` | genome 名称，例如 `hg38` |
| `chroms` | 染色体名 → 长度，顺序即文件顺序 |
| `bp_resolutions` | BP 分辨率 |
| `frag_resolutions`、`sites` | FRAG 分辨率及每条染色体排好序的酶切位点；有 FRAG 分辨率时必须提供 `sites` |
| `version` | 默认 9，也可写 8 |
| `attributes` | 键值属性，默认 `{"software": "omicsweave"}` |

两种 writer 共同的选项：`level` 是 zlib 压缩等级（0–9，默认 6，与 Juicer 工具相同）；`blocks=HicBlocks(intra, inter)` 设置块边长（默认染色体内 128 bin、染色体间使用 Juicer 网格；`"juicer"` 表示 Juicer 工具的网格，约 1000 bin）；`threads` 大于 1 时多线程编码和压缩块，输出与线程数无关；`sidecar=True` 同时写 [sidecar](sidecar.md)。

块大小只影响读取效率，任何读取方都从文件中获取块大小。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 提示分辨率不存在 | 用 `reader.resolutions()` 选择文件中已有的分辨率 |
| 矩阵全为 0 | 检查染色体名和 `reader.maps`；该染色体对可能没有接触 |
| 数值与 Juicebox 显示不同 | Juicebox 默认可能显示归一化值；OmicsWeave 返回原始计数 |
| 需要 KR / ICE 归一化或 O/E | 在返回的原始矩阵上自行计算，或使用专门的 Hi-C 分析工具 |
| pixels 的 row / col 很大 | 它们是 bin 起点的基因组坐标，减去窗口起点再除以分辨率得到下标 |
| 对角窗口出现重复像素 | 对称格子 `(i, j)` 与 `(j, i)` 都会返回 |

## 相关格式

- [cool](cool.md)、[mcool](mcool.md)：基于 HDF5 的接触矩阵，转换为 `.hic` 后查询。
- [4DN pairs](pairs.md)、[Juicer short](juicer-short.md)：每行一条接触的文本，转换时分箱。
- [Sidecar](sidecar.md)：记录 `.hic` 每块解压后的大小和记录数，供数据加载器预先规划，普通查询不需要。

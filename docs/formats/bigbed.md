# bigBed

[格式总览](README.md) · [查询](../query.md) · [转换](../convert.md) · [定量](../quant.md) · [BED](bed.md) · [bigWig](bigwig.md)

bigBed 是 [BED](bed.md) 的**带索引二进制形式**：它保存同样的区间和附加列，但把记录压缩成块、加上 R 树索引和描述各列的 schema，因此可以只读取与某个区域重叠的记录，而不必扫描整份文件。peaks、注释特征、甲基化调用、重复元件、黑名单区域等“区间 + 属性”数据，都适合存为 bigBed。

在 OmicsWeave 中，bigBed 是区间记录的标准查询格式：BED、bedGraph、bedMethyl 转换为 bigBed 后，由 `IntervalReader` 按区域读取；`quant` 的 `Peaks` 特征也直接读取 bigBed，把每条区间作为一个计数特征。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 基因组区间，及 BED 标准列和自定义附加列 |
| 常见后缀 | `.bb`、`.bigbed`、`.bigBed` |
| 存储 | 二进制 BBI 容器，数据块用 zlib 压缩 |
| 坐标 | 0-based 半开，与 BED 相同 |
| 随机访问 | 支持，文件内置 R 树索引 |
| 外部索引 | 不需要 |
| schema | header 记录字段数、标准字段数，可附带 autoSql |
| 典型来源 | UCSC `bedToBigBed`、ENCODE、UCSC 轨道下载、`ow convert` |
| OmicsWeave | 区域查询、直接写入、作为转换目标、作为 quant 特征 |

## 格式由来

bigBed 与 [bigWig](bigwig.md) 由 UCSC 基因组浏览器团队同时设计（Kent et al., *Bioinformatics* 2010），二者共用 **BBI**（Big Binary Indexed）容器：相同的 header、染色体树、R 树索引和 zoom 层结构，区别只在数据块里放的是数值信号还是区间记录。设计初衷是让浏览器只取可见区域，这与训练数据加载“给定区域取记录”的模式一致。

文本 BED 可以有任意自定义列，但没有地方说明这些列是什么。bigBed 用 header 中的字段计数和可选的 **autoSql** 描述补上了这一点。

## 文件结构

```text
┌──────────────────────────────┐
│ 固定 header                  │  magic、版本、字段数、标准字段数、各段偏移
│ zoom header 表及预留空间     │  各汇总层的分辨率与位置
│ autoSql 文本（可选）         │  每一列的名称、类型与说明
│ 全文件 summary               │
│ extension header             │  扩展头；OmicsWeave 不写 extra index
│ 染色体 B+ 树                 │  染色体名 → 编号、长度
│ 记录数（uint64）             │  原始 bigBed 记录总数
│ 数据块（zlib 压缩）          │  按染色体、起点排序的记录
│ R 树索引                     │  区域 → 数据块位置
│ zoom 层数据 + 索引           │  区间覆盖度的预汇总
│ 尾部 magic                   │
└──────────────────────────────┘
```

**数据块**中每条记录保存为：染色体编号、起点、右端，以及**第 4 列起的所有字段连成的一个 tab 分隔字符串**。也就是说，bigBed 并不把 name、score、strand 拆成二进制字段，而是原样保留文本。

**header 中的两个字段数**描述记录的结构：

- `fieldCount`：每条记录的总列数（含前三列坐标）。
- `definedFieldCount`：前多少列遵循 BED 标准含义，其余为自定义列。例如 narrowPeak 是 BED6+4：`fieldCount=10`、`definedFieldCount=6`。

**autoSql** 是 UCSC 用来描述表结构的简短语言，例如 bedMethyl 转换写入的 schema 开头为：

```text
table bedMethyl
"modkit pileup: modified base calls per reference position"
    (
    string chrom;             "Reference sequence name"
    uint   chromStart;        "Start position (0-based)"
    uint   chromEnd;          "End position (exclusive)"
    string name;              "Modified base code, with the motif when there are several"
    ...
    float  percent_modified;  "Nmod / Nvalid_cov * 100"
    ...
    )
```

它让下游工具知道每一列的名称和类型，但只是描述，记录本身仍以文本保存。

上图是 OmicsWeave writer 的段顺序；读取其他文件时应使用 header 中的偏移，不能假定这些段紧邻排列。

**zoom 层**保存区间覆盖度的预汇总，主要供浏览器在缩小视图时显示密度；OmicsWeave 默认按记录密度推导并保留能减少记录数的层，空文件及部分小文件没有 zoom 层。查询区间记录时不使用 zoom。

## 关键概念

### 返回完整记录

区间查询返回满足 `record.end > region.start` 且 `record.start < region.end` 的记录，空查询区域返回空结果。`fetch(region, clip=False)` 是默认行为：跨越查询边界的记录保留原始起止坐标。需要把坐标截到查询区域时传 `clip=True`。bigWig 的 `intervals` 默认则为 `clip=True`。

裁剪只改变返回的 `start`、`end`，附加字段保持原样，包括 BED12 的 `thickStart`、`thickEnd`、`blockSizes`、`chromStarts`。`clip=True` 不会重建一个裁剪后的 BED12 转录本结构。

### 区间可以重叠和嵌套

bigBed 允许任意重叠，包括长区间完全包住短区间。索引能正确找到这种长区间，不会因为它起点很早而漏掉。

### BED12 是整体记录

BED12 的 blocks（外显子结构）保存在附加字段中，索引只认识记录的整体跨度 `[chromStart, chromEnd)`。因此查询落在内含子里的区域，也会返回这个转录本。需要按外显子判断时，自行解析 blockSizes / blockStarts，或使用 [GTF](gtf.md) 的 `Annotation`。

### 链不参与查询

即使记录有 strand 列，按区域查询也不做链过滤；负链 key 返回的记录与正链相同，并保持基因组顺序。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 格式识别 | `formats.detect` 按 magic 识别为 `BigBed` |
| 区域查询 | `IntervalReader.fetch`：重叠记录的 `start`、`end`（`uint32`）与附加字段（`StringVector`），按起点排列 |
| schema | `IntervalReader.schema`：`field_count`、`defined_field_count`、`auto_sql`；`fields` 为解析后的字段，`columns(hits)` 按 autoSql 类型拆成各列数组 |
| 文件信息 | `chroms`：染色体名 → 长度的只读 mappingproxy |
| 缓存与并发 | 可接 `BlockCache`；解码释放 GIL；可 `clone()`、pickle |
| 写入 | `BigBedWriter`，可多线程压缩块，可选同时写 sidecar |
| 转换生成 | BED、bedGraph、bedMethyl → bigBed；MAF → bigMaf（也是 bigBed） |
| 定量 | `quant.Peaks("peaks.bb")`：每条记录一个特征列 |

限制：

- 查询 key 只接受连续 `Region`，不接受 `Spliced`。
- 不提供按名称查找（bigBed 的 extra index）和 zoom 层查询。
- 只读本地文件，只支持 little-endian 文件。
- 当前不提供 bigBed → BED 的文本导出。

## 查询

```python
from omicsweave.query import IntervalReader

with IntervalReader("peaks.bb") as peaks:
    print(dict(peaks.chroms))           # {'chr1': 248956422, ...}
    print(peaks.schema)                 # BigBedSchema(field_count=6, defined_field_count=6, auto_sql=None)

    hits = peaks.fetch(("chr1", 100, 200))
    print(len(hits))
    print(hits.start, hits.end)         # uint32 数组，0-based 半开，默认不裁剪
    print(hits.annotation.to_list())    # ['peak1\t960\t+', ...]
```

`annotation` 是一个 `StringVector`：`to_list()` 得到 Python 字符串列表，`bytes(i)` 取单条的原始字节，`to_arrow()` 转为 Arrow 字符串数组（需要 `pyarrow`）。需要某一列时按 tab 拆分：

```python
names = [row.split("\t")[0] for row in hits.annotation.to_list()]
scores = [int(row.split("\t")[1]) for row in hits.annotation.to_list()]

# 或按 schema 一次拆列：数值字段得到对应 dtype 的数组
cols = peaks.schema.columns(hits, ["name", "score"])
cols["score"]             # uint32 [N]
cols["name"]              # StringVector
cols["name"].to_list()    # list[str]
```

`schema.fields` 返回全部列的 `AutoSqlField(name, type, comment)`，前三列是坐标；`schema.columns(...)` 只拆分第 4 列起的注释字段，既接受 `hits`，也接受 `hits.annotation`。省略 `names` 时拆分全部附加列。

`AutoSqlField.dtype` 给出标量数值类型的 NumPy dtype，例如 `uint` → `uint32`、`int` → `int32`、`float` → `float32`、`double` → `float64`。字符串、数组、enum、set，以及颜色字段 `reserved` / `itemRgb` 的 dtype 为 `None`，拆列后得到 `StringVector`。拆列会分配新结果；行的字段数不符或数字无法转换时会报错。

没有 autoSql，或解析出的列数与 `field_count` 不符时，`fields` 回退到前 12 列的标准 BED 名称，其余列叫 `field13`、`field14` 等。此回退只是列名约定，不推断自定义列的真实语义；使用这类文件前应核对数据约定。

### 生成标签

区间数据常用来为训练样本生成标签，例如窗口内的 peak 覆盖：

```python
import numpy as np

def peak_mask(peaks: IntervalReader, chrom: str, start: int, end: int) -> np.ndarray:
    hits = peaks.fetch((chrom, start, end), clip=True)
    mask = np.zeros(end - start, dtype=np.uint8)
    for s, e in zip(hits.start, hits.end):
        mask[int(s) - start : int(e) - start] = 1
    return mask
```

`clip=True` 让坐标落在窗口内，可以直接换算成窗口下标。

## 准备 bigBed

### 从文本转换

```sh
ow convert peaks.bed peaks.bb --chrom-sizes genome.fa.fai
ow convert signal.bedGraph signal.bb --chrom-sizes genome.chrom.sizes
ow convert methylation.bed methylation.bb
```

转换会排序输入、允许重叠，并保留第 4 列起的所有字段。BED 和 bedGraph 转换不写 autoSql，标准字段数取 `min(列数, 12)`；可识别的 [bedMethyl](bedmethyl.md) 会保留全部 18 / 11 列，标准字段数为 9，并写入相应方言的 autoSql。各来源的规则见 [BED](bed.md)、[bedGraph](bedgraph.md)、[bedMethyl](bedmethyl.md)。

### 从内存写入

```python
from omicsweave.formats import BigBedSchema, BigBedWriter

chroms = {"chr1": 1000, "chr2": 800}
with BigBedWriter("example.bb", chroms, BigBedSchema(field_count=4), sidecar=True) as writer:
    writer.add("chr1", [100, 150], [200, 250], ["peak1", "peak2"])  # 每条一个字段串
    writer.add("chr2", 300, 400, "peak3")                           # 标量形式
```

- `chroms` 的顺序就是文件中的染色体顺序；记录必须按这个顺序、按起点排序。区间可以重叠。
- `annotation` 是第 4 列起的字段，多列用 `\t` 连接；可以给一个字符串（所有记录共用）、每条一个字符串的序列，或 `StringVector`。省略时为空。字段数应与 schema 一致。
- `BigBedSchema(field_count, defined_field_count=None, auto_sql=None)`：`field_count` 含前三列，默认 3；`defined_field_count` 默认等于 `field_count`。

带自定义列和 autoSql 的例子：

```python
auto_sql = """table scoredPeak
"Peaks with a signal value"
    (
    string chrom;      "Reference sequence"
    uint   chromStart; "Start (0-based)"
    uint   chromEnd;   "End (exclusive)"
    string name;       "Peak name"
    float  signal;     "Signal value"
    )
"""
schema = BigBedSchema(field_count=5, defined_field_count=4, auto_sql=auto_sql)
with BigBedWriter("scored.bb", {"chr1": 1000}, schema) as writer:
    writer.add("chr1", [100, 300], [200, 400], ["p1\t5.2", "p2\t1.7"])
```

`options=BbiOptions(...)` 可调整块大小（`items_per_slot`，bigBed 默认 512）、索引扇出、压缩等级和 zoom 层；`threads` 大于 1 时在工作线程上压缩块，输出与线程数无关。输入尚未排序时，可以先写成 BED，再用 `ow convert` 转换。

## 用于定量

```python
from omicsweave import quant

product = quant.Counts.chip(quant.Peaks("peaks.bb"))
```

```sh
ow quant atac.bam -o atac --atac --peaks peaks.bb
ow quant chip.bam -o chip --chip-peaks peaks.bb
```

`Peaks` 把每条 bigBed 记录作为一个特征列，特征 id 为 `chrom:start-end`。记录的完整区间参与计数，附加列中名称相同的记录不会被合并。BED 需先转换为 bigBed。用法见 [quant](../quant.md)。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 结果坐标超出查询区域 | bigBed 默认返回完整记录；需要时传 `clip=True` |
| 查询内含子也返回转录本 | BED12 按整体跨度索引；自行解析 blocks |
| `annotation` 只有一个字符串 | 第 4 列起的字段以 tab 连接；用 `schema.columns(hits)` 按列拆分 |
| `auto_sql` 为 `None` | 文件没有 autoSql（BED 转换不写），字段含义由 `defined_field_count` 和数据约定决定 |
| 找不到染色体 | 检查 `peaks.chroms`；转换时只写入有记录的染色体 |
| 写入时报排序错误 | `BigBedWriter` 要求按染色体顺序和起点排序；先排序，或改用 BED + `ow convert` |

## 相关格式

- [BED](bed.md)：bigBed 的文本来源。
- [bigWig](bigwig.md)：同一 BBI 容器中的数值信号。
- [bedMethyl](bedmethyl.md)：转为 bigBed 时保留全部列和 autoSql。
- [bigMaf](bigmaf.md)：bed3 + 整个 MAF 比对块的 bigBed，由 `MsaReader` 投影为比对矩阵。
- [Sidecar](sidecar.md)：bigBed 数据块布局的可选伴随文件，普通查询不需要。

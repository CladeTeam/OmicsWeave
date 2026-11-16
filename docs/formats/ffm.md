# FFM

[格式总览](README.md) · [转换](../convert.md)

> 状态：格式版本 1.0。Rust 读写（`crates/ffm`）、外部格式导入（`ow convert`，见 [Motif 文本格式](motif-texts.md)）与 Python 读写（`query.MotifReader`、`formats.FfmWriter`，值为 `types.Motif`）已实现；匹配见 [Motif 匹配](../motif.md)。

FFM（Factor Frequency Model）是 OmicsWeave 设计的 **motif 集合二进制格式**。一个 motif 被描述为若干「局部因子」：每个因子覆盖 motif 窗口内的少数几个位置，配一张非负频率表。PWM 只是每个因子覆盖一个位置的特例；相邻二核苷酸模型、非齐次 Markov 模型（BaMM）、带远距离位置对的势模型、k-mer 势表都能用同一结构表示。文件对已存储的 `f64` 无损；从对数权重转换成势值需要单独考虑浮点精度和范围。

所有元信息（字母表、注释键、ID、块目录、逐 motif 索引）集中在文件头部，打开时两次读取即可全部取得；块数据在其后按块压缩。按序号 `O(1)` 定位任意一个 motif，只需一次读取并解压它所在的块；不需要外部索引文件。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 一组 motif：每个 motif 由宽度、建模符号集合和若干因子频率表组成，另附可选注释 |
| 后缀 | `.ffm` |
| 存储 | 二进制，小端序；按块组织，块内分结构流、数值流、注释流三条流独立压缩 |
| 字母表 | 任意 1–255 个单字节符号，文件内统一定义；可选互补关系 |
| 数值 | IEEE 754 `f64`，有限且 ≥ 0 |
| 压缩 | 数值流 pcodec，注释流 zstd，也可不压缩 |
| 随机访问 | 按序号 `O(1)`、按 ID 期望 `O(1)` |
| 外部索引 | 不需要 |
| 不包含 | 背景分布、伪计数、阈值、p 值表、扫描策略 —— 这些是匹配层的参数 |

## 设计原则

- **只描述模型本身。** 文件里只有字母表、因子结构和频率表；外部格式里的 `E=`、`URL`、各数据库特有的表头都不进入模型区。外部格式只是 `ow convert` 的输入。
- **模型与注释分离。** 注释与模型同属一条 motif 记录，但存放在独立的注释流里：只读模型时注释流完全不解压。
- **低复杂度优先，数据结构保持干净。** 定位任何 motif 都不需要搜索或扫描：逐 motif 的定长索引在写入时生成。能从元信息以 `O(n)` 廉价推导的加速结构（如 ID 哈希表）不写进文件，由读取方打开时在内存中构建。
- **存储忠实，解释后置。** 文件保存频率（可以是未归一化的计数），不预先换算成打分：伪计数、背景、整数化都是扫描时的选择，预先换算会丢信息。
- **紧凑的符号表示。** 每个 motif 只为它关心的符号建表，其余符号合并成一类，蛋白等大字母表下表的规模大幅缩小。

## 概念

### 字母表、编码与 invalid

字母表是文件级的：1–255 个**符号**，每个符号是一个字节，按顺序编号为**全局编码** `0 … n-1`。字母表还给出一张 256 项的**查找表（LUT）**，把任意输入字节映射到全局编码；同一编码可以有多个字节（例如 `a` 与 `A`，或把 `U` 映射到 `T`）。

LUT 中值为 `0xFF` 的字节是 **invalid**：它不属于字母表，代表「未知」而不是某个符号。例如 DNA 字母表 `ACGT` 下的 `N`、IUPAC 模糊碱基和 gap 都应映射为 invalid。含 invalid 的窗口如何处理（跳过、报错等）由扫描时的策略决定，格式不规定。

可选的**互补表**为每个全局编码给出它的互补编码，必须是对合（`comp[comp[c]] == c`）。有互补表的字母表才能扫描反向链。

### support 与 rest

每个 motif 有自己的 **support**：它单独建模的全局编码集合（`s` 个，`1 ≤ s ≤ n`）。字母表中不在 support 里的符号合并成一类 **rest**。

| 情形 | 局部编码 | 局部基数 `r` |
|---|---|---|
| `s = n`（DNA motif 的常态） | support 按全局编码升序编为 `0 … n-1`，没有 rest | `r = n` |
| `s < n` | support 按全局编码升序编为 `0 … s-1`，rest 编为 `s` | `r = s + 1` |

rest 是一个真实的类别，在频率表中占一格：它的频率是「这个位置出现 support 以外任何符号」的总频率。invalid 不属于 rest，也不在任何表中出现。

例：蛋白字母表（`n = 20`）上一个只关心 `D E K R` 的 motif，`s = 4`、`r = 5`，一个二元因子的表只有 25 格，而按全字母表需要 400 格。

### 因子

motif 有宽度 `W`（位置 `0 … W-1`）。一个**因子**由三部分组成：

- **类型**：`joint`、`conditional` 或 `potential`；
- **位置**：`k` 个严格递增的位置 `p_1 < … < p_k`（`1 ≤ k ≤ 16`），称 `k` 为因子的阶；
- **表**：`r^k` 个 `f64`，下标为局部编码按位置顺序排列、最后一个位置变化最快：

  `index = Σ_j code(p_j) · r^(k-j)`

| 类型 | 表的含义 | 归一化组（和必须 > 0） |
|---|---|---|
| `joint` | 这些位置上局部编码组合的联合频率 | 整张表 |
| `conditional` | 给定前 `k-1` 个位置时，最后一个位置的条件频率 | 每个上下文对应的连续 `r` 项 |
| `potential` | 未归一化的势：打分贡献为 `ln(值)` | 无；值为 0 表示禁止该组合 |

`joint` 与 `conditional` 的表可以是未归一化的计数，读取方按归一化组相除得到频率。没有被任何因子覆盖的位置不贡献信息；同一位置可以被多个因子覆盖。

### sites

motif 可选地携带 `sites`（`f64`，有限且 > 0）：频率背后的观测数，例如 JASPAR 计数矩阵的列和或 MEME 的 `nsites`。伪计数平滑需要它。没有 `sites` 时，由匹配层采用默认值。

### 带 gap 的 profile

motif 可选地携带 gap 计数，成为**带 gap 的 profile**（GLAM2 一类的比对模型）：`W` 个位置是比对的列，序列可以在列之间插入残基，也可以缺失某一列。此时因子必须恰好是每列一个单位置 `joint` 因子、按列顺序排列，表为该列的残基计数；另有两组计数：

- `deletions`：`W` 个值，每列被缺失的序列数；
- `insertions`：`W − 1` 个值，第 `i` 列与第 `i+1` 列之间插入的残基总数。

计数有限且 ≥ 0；`sites` 通常是比对的序列数。`deletions` 大于 0 的列允许残基计数全为 0（所有序列都缺失该列），这是对「归一化组之和必须 > 0」的唯一例外。带 gap 的 profile 由比对打分，不按因子模型扫描。

### 常见模型的表示

| 来源 | 因子 | `sites` |
|---|---|---|
| PFM / 计数矩阵（JASPAR、TRANSFAC、HOCOMOCO PCM） | `W` 个一阶 `joint`，表为计数 | 列和 |
| 概率矩阵（MEME、HOMER、CIS-BP） | `W` 个一阶 `joint`，表为概率 | MEME 的 `nsites`；没有则缺省 |
| 相邻二核苷酸模型（HOCOMOCO di 模型） | `W-1` 个二阶 `joint`，位置 `(i, i+1)` | 有则保留 |
| k 阶非齐次 Markov（BaMM） | 每个位置 `i` 一个 `conditional`，位置 `max(0, i-k) … i` | 有则保留 |
| 能量模型（ProBound） | 单位置与位置对的 `potential`，值为 `exp(权重)` | — |
| k-mer 打分表（PBM 等） | 一个 `k` 阶 `potential`，位置 `0 … k-1` | — |
| GLAM2 比对 | 每列一个一阶 `joint`（残基计数）+ gap 计数 | 比对序列数 |

## 文件结构

```text
┌────────────────────────────────────────────┐
│ Header（128 字节，定长）                     │ magic、版本、各段位置
├────────────────────────────────────────────┤  ┐
│ Alphabet                                    │  │
│ Schema                                      │  │ 元信息：打开时
│ IDs                                         │  │ 一次读取
│ Directory                                   │  │
│ Index                                       │  │
├────────────────────────────────────────────┤  ┘
│ Block 0                                     │ 结构流 + 数值流 + 注释流
│ Block 1                                     │
│ …                                           │
└────────────────────────────────────────────┘
```

所有整数与浮点数均为小端序。校验使用 `xxh3_64`：

- Alphabet、Schema、Directory 以 8 字节校验和结尾，覆盖该段内校验和之前的所有字节；
- 每条**流**（IDs、Index 和块内的三条流）自带校验和，见下文「流描述」。只解码一条流时只校验这条流。

元信息依次相接，紧跟 Header；块数据从 Index 之后开始，直到文件末尾。写入方在所有 motif 写完后才能确定 ID、目录与索引，因此先在内存中保留各块的压缩结果，`finish` 时依次写出 Header、元信息和块。块是压缩后的字节，所需内存约等于最终文件大小。

### 复杂度

| 操作 | 代价 |
|---|---|
| 打开 | 2 次读取（Header；整段元信息）；解码、校验并构建 ID 哈希表，代价与元信息字节数、motif 数和块数相关 |
| 按序号取 motif | `O(1)` 定位：Index 第 `i` 项给出块号与流内偏移；块未缓存时一次读取整块并解码 |
| 按 ID 取 motif | 期望 `O(1)`：查读取方构建的哈希表，再按序号取 |
| 取注释 | `O(1)` 定位；首次访问需解码并校验整个块的注释流，数值字段解析另与文本长度相关 |
| 全量遍历 | 每块读取、解码各一次 |

文件读取接口 `FfmReader::motif(i)` / `annotation(i)` 返回 `Result`，序号越界返回 `Err`；下面常驻集合 `MotifSet` 的同名方法返回 `Option`，越界为 `None`。

### 常驻内存接口

motif 文件较小时，推荐 `MotifSet::load(source)`：加载时读取所有块、完成模型和注释校验，返回与输入源独立的只读集合。`FfmReader::into_set()` 可将已打开的 reader 消费为同样的集合。

```rust
use omicsweave_core::FileSource;
use omicsweave_ffm::MotifSet;

let motifs = MotifSet::load(FileSource::open("motifs.ffm")?)?;
if let Some(entry) = motifs.get_by_id("MA0001.1") {
    let model = entry.model();
    let factor = model.factor(0).unwrap();
    let table: &[f64] = factor.table;
    let name = entry.annotation().get_by_name("name");
    if let Some(matrix) = model.pfm() {
        let first_position: &[f64] = matrix.row(0).unwrap();
    }
}
```

| 对象 / 接口 | 说明 |
|---|---|
| `MotifSet` | 共享所有权的常驻集合；`clone()` 共享内存，查询不加锁、不分配、不访问文件 |
| `get(i)` / `get_by_id(id)` / `iter()` | 借用的 `MotifEntry`，同时访问 ID、模型、注释 |
| `motif(i)` / `annotation(i)` | 单独访问模型或注释；越界返回 `None` |
| `MotifView::factor(i)` | 常驻集合与文件块 PFM 视图内 `O(1)` 获取因子；其他文件块视图内为 `O(i)` |
| `table_values()` | 所有表连续存放的 `&[f64]`；不含 `sites`；集合与单 motif 均可调用 |
| `MotifView::pfm()` | 恰好每个位置一个一阶 joint 因子、按位置排列时，返回零拷贝 `PfmView` |
| `PfmView::values()` / `row(p)` / `rows()` | 按位置优先的连续矩阵，列为局部编码，包括存在时的 rest |
| `AnnotationView::get(i)` / `value(i)` | 字段偏移直接索引；后者按 schema 解释为 `AnnotationValue`，缺失返回 `None` |
| `heap_bytes()` | 当前集合持有的堆内存，包括数组预留容量，不含共享引用计数与 allocator 管理开销；共享 clone 不重复分配 |
| `write_to(out, options)` | 直接从借用视图写文件，不构造嵌套 `Vec<Factor>` |

集合内部使用 `MotifArena`：motif 描述、一般因子的描述、位置、support、`f64` 表、注释文本与字段偏移各自集中存放。PFM 直接表示为连续矩阵，不为每个位置重复保存因子描述；所有 PFM 共享一张 `0 … 最大宽度-1` 的位置数组。加载时建立偏移并选定布局，查询时无需重读结构流，冻结集合时收紧数组的预留容量。`f64` 数组满足 Rust 的自然对齐要求，不承诺额外的 SIMD 对齐。注释在常驻集合中全部解码；`FfmReader` 仍按需解码注释。

视图只借用集合，不能比所借用的集合句柄活得更久。需要跨线程共享整个集合时 clone `MotifSet`；需要独立编辑一个模型时使用 `to_motif()`，此时才为各因子分配独立数组。原有文件读取接口的 `MotifRef` / `AnnotationRef` 持有块的共享所有权，可在 reader 释放后继续使用。

程序中直接构造集合可使用 `MotifSetBuilder::push(id, &model, annotation)`，或 `push_pfm(id, values, sites, annotation)` 从扁平的全字母表 PFM 构造，避免先建立每行一个 `Vec`。`push_view` 可直接复制已校验的借用视图；其全局编码必须与目标字母表使用相同顺序。

### Python 读写与所有权

高层接口为 `formats.FfmWriter(path, alphabet, schema=None, *, block_records=4096, block_values=4096)` 与 `query.MotifReader(source)`：

```python
from omicsweave.formats import FfmWriter
from omicsweave.query import MotifReader
from omicsweave.types import Alphabet, Motif

dna = Alphabet.dna()
model = Motif.from_pfm("AC", [[9, 1, 1, 1], [1, 9, 1, 1]], annotation={"source": "demo"})
with FfmWriter("models.ffm", dna, {"source": "string"}) as writer:
    writer.add(model)

with MotifReader("models.ffm") as reader:
    ids = reader.ids                  # tuple[str, ...]，按文件顺序
    model = reader["AC"]              # 也支持序号、负序号、slice、迭代
    matrix = model.matrix             # 新的 (width, alphabet_len) float64 矩阵
    width = reader.widths()            # uint32[motif_count]
# 已取出的 model / matrix 数据在 reader 关闭后仍然有效。
```

`schema` 是保持顺序的 `name -> "string" | "integer" | "float"` 映射；`add` 接受 `Motif`，`add_all` 接受其迭代器。缺失注释写为空字段，读取为 annotation 中没有该键；原始字符串、整数、浮点按 schema 转为 Python 值。写入器关闭时输出文件；关闭前须保留完整的压缩块集合。`block_values` 是分块目标值，单个 motif 不跨块，因此一条大记录可以超过该目标。

Python `MotifReader` **整体载入到 Rust `MotifSet`**，查询不访问文件；它没有文件块 `BlockCache` 参数。`reader[i]` 返回的 `types.Motif` 持有复制到 NumPy 的因子表，不是 Rust 借用视图；反复取同一模型会再次分配。编译整个集合时直接将 reader 传给 `motif.compile`，可避免逐个构造 Python 模型。Rust `FfmReader` 的按需解码/显式 `BlockCache` 与 Python 常驻 reader 是两种不同接口。

`Factor.table` 是含 `arity ** order` 格的一维 `float64` 数组，`Gaps.deletions/insertions` 分别是长 `width` / `width-1` 的一维 `float64` 数组；positions 是递增的 motif 位置，表最后一个位置的编码变化最快。`support` 为全局编码的 bytes 或 `None`；存在子集时表的 arity 包含最后一个 rest 列。`Motif.matrix` 只为全字母表、每个位置一个 joint 因子的模型返回矩阵；一般因子模型或 support 子集返回 `None`。构造 Python 值时 `ascontiguousarray` 可以复用调用方数组，frozen dataclass 不等于底层数组只读；写入与编译会在原生边界校验模型。

读取接口见 [`query/_motif.py`](../../src/omicsweave/query/_motif.py)，写入接口见 [`formats/_motif.py`](../../src/omicsweave/formats/_motif.py)，值对象见 [`types/_motif.py`](../../src/omicsweave/types/_motif.py)。

### 读取资源限制

`FfmReader::open_with_options` 和 `MotifSet::load_with_options` 接受 `ReaderOptions`。默认元信息、单块压缩字节和每条流的解码字节分别最多 256 MiB，常驻逻辑载荷最多 1 GiB。高层 Python reader 使用这些默认值，当前不暴露 `ReaderOptions` 配置参数。它们是可调整的读取策略，不改变 1.0 文件布局。

读取方在分配前检查长度与限制；索引流必须先满足 `decoded_len = 4 × motif_count`。常驻限制覆盖模型描述、表、位置、注释偏移及元信息；输入源、临时解码缓冲和数组预留容量还会增加峰值内存，不能把该限制当作进程 RSS 上限。

### Header

| 偏移 | 长度 | 字段 | 说明 |
|---:|---:|---|---|
| 0 | 8 | magic | `89 46 46 4D 0D 0A 1A 0A`（`\x89FFM\r\n\x1a\n`） |
| 8 | 2 | major | `1`；读取方拒绝不认识的 major |
| 10 | 2 | minor | `0`；minor 增加只引入向后兼容的扩展 |
| 12 | 4 | flags | 保留，必须为 0；读取方拒绝不认识的位 |
| 16 | 8 | motif_count | motif 总数，受 Index 和 IDs 的流长度限制，见「限制」 |
| 24 | 8 | block_count | 块数 |
| 32 | 16 | alphabet | 段的偏移（固定为 128）与长度，各 `u64` |
| 48 | 16 | schema | 偏移与长度 |
| 64 | 16 | ids | 偏移与长度 |
| 80 | 16 | directory | 偏移与长度 |
| 96 | 16 | index | 偏移与长度 |
| 112 | 8 | reserved | 必须为 0 |
| 120 | 8 | checksum | `xxh3_64(bytes[0..120])` |

各段长度都包含其校验和或流描述。五个元信息段按上表顺序依次相接；块数据区从 Index 末尾开始，到文件末尾结束。magic 的首字节不是 ASCII，`\r\n` 与 `\x1a` 能暴露文本模式传输造成的损坏，与 PNG 的做法相同。

### Alphabet

| 字段 | 类型 | 说明 |
|---|---|---|
| symbol_count | `u8` | `n`，1–255 |
| flags | `u8` | bit 0：有互补表；其余位必须为 0 |
| name_len | `u16` | |
| name | `name_len` 字节 | UTF-8，可为空；仅供显示，例如 `DNA` |
| symbols | `n` 字节 | 每个全局编码的规范字节 |
| complement | `n` 字节 | 仅当 flags bit 0 置位；每个编码的互补编码 |
| lut | 256 字节 | 字节 → 全局编码，`0xFF` 为 invalid |
| checksum | `u64` | |

约束：`symbols` 互不相同；`lut[symbols[c]] == c`；`lut` 的值只能是 `0 … n-1` 或 `0xFF`；互补表是对合。

OmicsWeave 会提供几个预置字母表（DNA、RNA、蛋白），但它们只是 API 的便利，文件里总是写出完整定义。

### Schema

| 字段 | 类型 | 说明 |
|---|---|---|
| key_count | `u16` | 注释键的个数，可以为 0 |
| 每个键：type | `u8` | `0` 字符串、`1` 整数（十进制 `i64`）、`2` 浮点（十进制 `f64`） |
| 每个键：name_len | `u16` | |
| 每个键：name | `name_len` 字节 | UTF-8，非空，不含 `\t` `\n` `\r` `\0`，互不相同 |
| checksum | `u64` | |

motif ID 不是注释键，单独存放在 IDs 段。

### 流描述

IDs、Index 两段与块内的三条流都以同样的 32 字节描述开头，紧跟流的数据：

| 偏移 | 字段 | 类型 | 说明 |
|---:|---|---|---|
| 0 | codec | `u8` | `0` 不压缩、`1` pcodec（standalone 格式）、`2` zstd |
| 1 | reserved | 7 字节 | 0 |
| 8 | decoded_len | `u64` | 解码后的元素数（`u32` / `f64` 的个数，或字节数），`< 2^32` |
| 16 | stored_len | `u64` | 数据在文件中的字节数 |
| 24 | checksum | `u64` | `xxh3_64(描述的前 24 字节 ‖ 数据)` |

数值流（元素为 `u32` 或 `f64`）允许 codec `0` / `1`，字节流允许 `0` / `2`。

### Block

一个块容纳连续的若干条 motif 记录；**一条记录不跨块**。读取时一次读入整块，再在内存中分别解压各流：

```text
record_count   u32    ≥ 1，等于 Directory 中的值
reserved       u32    0
结构流         描述 + 数据（u32）
数值流         描述 + 数据（f64）
注释流         描述 + 数据（字节）
```

结构流与数值流在读入后立即解码；注释流的压缩字节留在内存中，第一次访问注释时才解压。

#### 结构流（`u32` 序列）

每条记录依次写出：

```text
width                 W ≥ 1
support_count         s，1 ≤ s ≤ n
support[s]            仅当 s < n：严格递增的全局编码
flags                 bit 0：有 sites；bit 1：有 gap 计数；其余位为 0
factor_count          F ≥ 1
每个因子：
  kind                0 joint、1 conditional、2 potential
  order               k，1 ≤ k ≤ 16
  positions[k]        严格递增，< W
```

#### 数值流（`f64` 序列）

每条记录依次写出：可选的 `sites`；有 gap 计数时接着写 `W` 个 `deletions` 和 `W − 1` 个 `insertions`；然后按因子顺序写出每张表的 `r^k` 个值。所有值必须有限且 ≥ 0，joint / conditional 的归一化组必须有正质量；带 gap 且 deletion > 0 的空残基列除外，potential 不要求正质量；`sites` 必须 > 0。

#### 注释流（字节）

每条记录写出 `key_count` 个值，以 `\t` 分隔、以 `\0` 结尾。值不得包含 `\t` 与 `\0`；空串表示缺失。整数与浮点键的非空值必须能解析为对应类型；float 注释只要求可解析为 `f64`，不套用模型数值的有限非负约束。`key_count = 0` 时注释流为空（`decoded_len = 0`）。

### Directory

`block_count` 个 16 字节的定长项，按文件位置升序，之后是段校验和：

| 字段 | 类型 | 说明 |
|---|---|---|
| offset | `u64` | 块在文件中的起始位置 |
| length | `u32` | 块的字节数 |
| record_count | `u32` | 块内记录数，≥ 1 |

块从数据区起点开始依次相接，直到文件末尾；`record_count` 之和等于 `motif_count`。

### Index

一条 `u32` 流，`decoded_len = 4 × motif_count`，按列存放：

```text
block[motif_count]          所在块号
structure_at[motif_count]   在该块结构流中的起始元素
values_at[motif_count]      在该块数值流中的起始元素
annotation_at[motif_count]  在该块注释流中的起始字节
```

motif 按序号依次落在各块中：`block` 单调不减，每块恰好 `record_count` 项；同一块内三个偏移严格递增（无注释时 `annotation_at` 全为 0），块内第一条记录的偏移为 0。一条记录在某条流中的终点是同块下一条记录的起点，块内最后一条则是该流的长度。按列存放让每一列都是单调或锯齿状的整数，pcodec 压缩效果好。

### IDs

一条字节流，解码后的内容：

```text
offsets[motif_count + 1]   u32   motif i 的 ID 为 bytes[offsets[i] .. offsets[i+1]]
bytes                            所有 ID 依次相接
```

ID 必须是非空 UTF-8，不含 `\t` `\n` `\r` `\0`，在文件内唯一；ID 字节总长 `< 2^32`。文件不存储查找结构：读取方打开时据此构建哈希表。

## 读取

- **打开**：读 Header，再一次读取整段元信息；校验各段，检查 Index 与 Directory 一致，并用 IDs 构建内存哈希表（同时检查 ID 唯一）。
- **按序号取**：Index 第 `i` 项给出块号与偏移 → 块未缓存时一次读取整块、解码结构流与数值流、校验块内全部记录 → 直接切出这条记录。
- **按 ID 取**：查哈希表得到序号，再按序号取。
- **取注释**：同一块的注释流在第一次访问时解压，之后按 `annotation_at` 直接切出。
- 通过 `set_cache(Some(Arc<BlockCache>))` 显式共享解码后的块，与其他 reader 一致。未挂载缓存时，重复随机调用会各自解码；`motifs()` 在一次顺序遍历内复用当前块。需要反复访问整套模型时可加载 `MotifSet`。

## 参考语义

格式只规定表的含义；把频率换算成打分属于匹配层。作为参考，给定 0 阶背景 `b`（rest 的背景为 rest 中各符号背景之和）与平滑后的频率 `q`：

| 类型 | 打分贡献 |
|---|---|
| `joint`，位置 `P` | `ln q(x_P) − Σ_{p∈P} ln b(x_p)` |
| `conditional` | `ln q(x_last ∣ x_context) − ln b(x_last)` |
| `potential` | `ln(值)` |

motif 在一个窗口上的打分是所有因子贡献之和。PWM 时即为常规的 log-odds 打分。

匹配层先从原模型估计并平滑因子，展开到全局字母表，再在这些因子上构造反向链：位置 `p` 映射为 `W-1-p`，因子的位置轴反转、符号按全局互补表置换；之后才分别生成两条链的打分表与零分布。FFM 的 `MotifView::encoding(alphabet)` 则提供正反向两张字节查找表：先按全局互补表转换，再映射到 motif 的局部编码，invalid 始终保持 `0xFF`；反向读取时还需反转位置顺序。

不能普遍假定原 support 下的 rest 与自身互补。例如 DNA support 为 `{A}` 时，T 和 C 都是 rest，但互补后 T→A 属于 support，C→G 仍是 rest。若要把反向模型表示为重新排列后的局部表，需要同时把 support 转换成其互补集合并重新编号。只有 support 对互补封闭时，才可沿用原 support 并让 rest 映射到 rest。`conditional` 表还涉及条件方向，因此不直接在存储的频率表上做反向归一化。

`FactorView::normalize(arity, out)` 提供 joint / conditional 的稳定归一化：先按组内最大值缩放再求和，有限计数的总和超过 `f64::MAX` 时仍能正确归一化。potential 不适用此操作。

## 限制

| 项目 | 上限 |
|---|---|
| 字母表符号数 | 255 |
| motif 数 | 至多 `floor((2^32 − 1) / 4)`；另受 IDs 整条流的字节上限约束 |
| motif 宽度 | `2^32 − 1`；稀疏位置可以有巨大跨度，但匹配层另有宽度与上下文限制 |
| 因子阶 | 16 |
| 单张表的格数 `r^k` | `2^24` |
| 注释键数 | 65535 |
| 单块字节数 | `2^32 − 1` |
| 单块每条流的元素数 | `2^32 − 1` |

IDs 的 offsets 与文本合起来也是一条字节流，因此 `(motif_count + 1) × 4 + ID 字节总长 ≤ 2^32 − 1`。写入方在接受记录前检查这一条件，以及块内结构、数值和注释流的各自长度；压缩失败不会提前占用新 ID，也不会丢弃已接受记录。

## 校验规则汇总

读取方遇到以下任一情况即报错，不做修补：magic、版本或 flags 不符；任一段或流的校验和不符；流解码后的元素数与 `decoded_len` 不符；字母表约束不满足；块未依次相接或未覆盖全部 motif；Index 与 Directory 或块内记录不一致；结构流或数值流长度与记录描述不一致；位置越界或不严格递增；模型因子或 gap 计数为负、NaN 或无穷；joint / conditional 归一化组没有正质量且不属于合法空残基列；`sites` ≤ 0；注释值个数不等于 `key_count`；ID 重复或含非法字符。

## 块大小、验证与性能测量

块越大，单次随机访问需要解码的数据通常越多；压缩比例也依赖实际表值与编码器，不能由块大小直接保证。默认 `WriterOptions` 使用 pcodec 编码数值/索引流、zstd 编码字节流；Rust 可选 raw 编码，Python 写入器当前只暴露分块参数。

优先保留来源的整数计数；不要为压缩而把概率任意取整。格式接受计数与概率，压缩收益取决于数值分布。有限非负计数即使普通 `f64` 求和溢出，也可通过 `FactorView::normalize` 稳定归一化；非法负值、NaN、无穷仍被拒绝。

[`crates/ffm/benches/access.rs`](../../crates/ffm/benches/access.rs) 提供合成集合的读写、显式缓存随机访问、常驻访问与堆载荷测量入口：`cargo bench -p omicsweave-ffm --bench access`。输入为内存中的字节，不能将其延迟解释成磁盘或远程 I/O 性能。比较时分别记录无缓存 reader、同一 caller-owned `BlockCache` 的热访问与 `MotifSet`；最终常驻载荷不等于加载过程的 RSS 峰值。

[`tests/integration.rs`](../../crates/ffm/tests/integration.rs) 注册 [`cases/serialization.rs`](../../crates/ffm/tests/cases/serialization.rs)、[`cases/resident.rs`](../../crates/ffm/tests/cases/resident.rs)、[`cases/corruption.rs`](../../crates/ffm/tests/cases/corruption.rs) 等场景，定义读写往返、文件/常驻模型对照、所有权、缓存与畸形输入的检查。运行结果需要当前 checkout 的命令和日志支持；本页不把历史测试数量、旧基准或无原始记录的外部数据数字当成当前结论。

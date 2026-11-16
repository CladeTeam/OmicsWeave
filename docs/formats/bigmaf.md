# bigMaf

[格式总览](README.md) · [查询](../query.md) · [转换](../convert.md) · [MAF](maf.md) · [bigBed](bigbed.md)

bigMaf 是 UCSC 为在浏览器中按区域显示多序列比对而定义的格式：它就是一个 [bigBed](bigbed.md)，每条记录对应 [MAF](maf.md) 中的一个比对块，区间是参考物种的坐标，附加的一列保存整个比对块的文本。因此它继承了 bigBed 的 R 树索引，可以只读取与查询区域重叠的块。

在 OmicsWeave 中，`ow convert` 从 MAF 生成 bigMaf；专用的 Rust `MsaReader` 读取它，把比对块投影到参考坐标，返回 `[物种, 长度]` 的 `uint8` 矩阵。Rust 层复用 bigBed 的索引和解码缓存，直接扫描解压后的字节完成投影；Python 负责参数检查和结果包装，解码与投影释放 GIL。`MsaReader` 要求 bigBed 有 4 列，且 autoSql 的第四个字段实际声明为 `string mafBlock` 或 `lstring mafBlock`；只在描述文字中出现 `mafBlock` 不算。普通 bigBed 不符合这个布局时不能作为 MSA 读取。`query.open(...)` 按 bigBed 类型返回 `IntervalReader`，需要矩阵时应显式创建 `MsaReader`。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 多序列比对块，按参考坐标索引 |
| 常见后缀 | `.bigmaf`、`.bigMaf`、`.bb`、`.bigbed` |
| 容器 | bigBed（`bed3+1`），autoSql 为 UCSC 的 `bigMaf.as` |
| 坐标 | 参考物种的 0-based 半开坐标 |
| 随机访问 | 支持，bigBed 内置索引 |
| 典型来源 | UCSC `mafToBigMaf` + `bedToBigBed`、`ow convert` |
| OmicsWeave | 作为转换目标；`MsaReader` 读为比对矩阵，`IntervalReader` 读原始记录 |

## 记录结构

autoSql：

```text
table bedMaf
"Bed3 with MAF block"
    (
    string chrom;      "Reference sequence chromosome or scaffold"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    lstring mafBlock;   "MAF block"
    )
```

每条记录：

| 列 | 内容 |
|---|---|
| `chrom` | 参考行 `src` 去掉物种前缀，如 `hg38.chr22` → `chr22` |
| `chromStart` / `chromEnd` | 参考行的 `start` 与 `start + size` |
| `mafBlock` | 整个块，每行以 `;` 结尾 |

`mafBlock` 按 kent `mafWriteDelimiter` 的写法重新生成，而不是照抄原文：

```text
a score=25976.000000;s hg38.chr22             10510000 2 +  50818468 GA;s panTro4.chrUn_GL393523    88624 2 +    405060 GA;i panTro4.chrUn_GL393523 N 0 C 0;...
```

- `a` 行只保留 `score`，固定六位小数；其他属性丢弃。
- `s`、`e` 行各字段按块内最宽的值对齐（名称左对齐，数字右对齐），`q`、`i` 行紧跟在对应 `s` 行之后。
- `src` 为 `.` 的行改名为 `defaultName`；`s` 行碱基数为 0 且带 `i` 行时写成 `e` 行（kent 的行为）。

`ow convert` 采用上述文本布局和 autoSql 定义，染色体长度取参考行的 `srcSize`，按染色体首次出现顺序和起点排序。支持的输入范围和检查见 [MAF](maf.md)。

## 读取为比对矩阵

```python
from omicsweave.query import MsaReader

species = ["hg38", "panTro4", "rheMac3", "mm10", "canFam3"]
with MsaReader("multiz100way.bigmaf", species, fill="-") as msa:
    aln = msa.fetch(("chr22", 10_600_000, 10_600_020))
    print(aln.species)        # 与 species 相同
    print(aln.bases.shape)    # (5, 20)
    for name in aln.species:
        print(name, aln.row(name).tobytes().decode())
```

投影规则：

| 情况 | 结果 |
|---|---|
| 参考行为 `-` 的列（其他物种的插入） | 丢弃，每列对应一个参考碱基 |
| 其他物种行为 `-` | 保留 `-`，表示相对参考的缺失 |
| 物种在块中没有 `s` 行（包括只有 `e` 行） | `fill` |
| 没有任何块覆盖的位置 | 所有行为 `fill` |
| 同一块中某物种有多行 | 取第一行 |
| 参考区间重叠的块 | 按记录顺序覆盖各块实际出现的物种行；后块缺失的物种不清除前块已有碱基 |
| key 在 `-` 链 | 列顺序反转，比对碱基按 IUPAC 互补表变换；`-` 与未覆盖位置的 `fill` 不变 |

`fill` 默认是 `"-"`，即“无比对”与缺失不区分；需要区分时可设为 `"."`。它可以是一个 ASCII 字符、长度为 1 的 bytes，或 0…255 的整数；bytes/整数也可指定非 ASCII 填充值，返回数组仍为原始 `uint8`。碱基保留文件中的大小写。负链查询中，未覆盖位置的填充值保持原样；例如 `fill="A"` 仍填 `A`，比对中实际出现的 `A` 则互补为 `T`。

物种名当前按 `src` 的第一个点切分。因此转换器虽能通过 `ref_species="hg19.1"` 识别含点的参考名，reader 中该行仍叫 `hg19`。不传 `species` 时，行为本次查询读到的物种，按首次出现排序；没有块的查询返回零行。固定 `species` 后，每次都返回给定的行数和顺序，缺失物种用 `fill` 补齐。

查询只接受连续 `Region`，不接受 `Spliced`。区间超出染色体末端时，仍返回完整的请求长度，末端以外的位置填 `fill`。`out=` 需要打开 reader 时固定 `species`，并传入可写、C-contiguous、形状为 `(len(species), len(region))` 的 `uint8` 数组；返回结果引用这个数组。省略 `out` 时分配独立结果。投影会逐块解析 MAF 文本，即使只选择少数物种也仍需检查读到的行；`cache=BlockCache(...)` 缓存的是底层 bigBed 解码块，不缓存投影后的矩阵。

### 校验、错误与所有权

输入以每个记录的第一个 `s` 行作为参考，要求它在 `+` 链；参考的 start/size、srcSize 与 bigBed 的完整区间、染色体长度一致，源名称须等于染色体名或以 `.染色体名` 结尾。`s` 行接受空白或 tab 分隔；所有行均检查七个字段、UTF-8 源名称、可打印 ASCII 文本、整数/链/范围、相同列数，以及 `size` 等于非 gap 字节数。固定少数物种只减少输出，不跳过其他 `s` 行的结构检查。参考负链记录不会被静默按正链投影，而是报错。

默认 `fetch` 分配独立、NumPy 自己拥有内存的二维数组（`OWNDATA=True`），后续查询、关闭 reader 或清空缓存不会覆盖已返回的数据。原生层在同一 reader 锁内完成 query 与 write，通过 NumPy C API 创建输出，不在持锁时调用可被 Python 替换的 `numpy.empty`。解析、投影和输出的大块写入释放 GIL；数组创建和包装需要 GIL。`out=` 则保留调用方数组身份，Python 与原生入口都检查完整二维形状，不只检查元素总数。

每个执行上下文复用 `bases`、列映射与覆盖来源缓冲。`clone()` 共享不可变 bigBed 元信息和已挂载的外部 `BlockCache`，保留相同缓存 namespace，但有独立的解码/投影 scratch；不重新读取 autoSql，也不复制上次查询状态。缓存 lease 保证 visitor 借用字节的生命周期，矩阵不借用解压块。Rust `query` 出错后 `write` 返回 `NotQueried`，下次成功查询可继续使用。

未知染色体与没有命中块相同：固定物种时返回请求长度的 fill，动态物种时返回 `(0, 请求长度)`。空区间的列数为 0；输出行列乘积必须适合地址空间。默认结果至少需要一份暂存矩阵和一份独立输出；fill 需要按来源区分互补时还会保留覆盖标记，不能把缓存上限当作整个查询的内存上限。

Rust 入口见 [`crates/msa/src/lib.rs`](../../crates/msa/src/lib.rs)，Python 绑定见 [`crates/python/src/msa.rs`](../../crates/python/src/msa.rs)；性质测试与绑定边界分别位于 [`projection_properties.rs`](../../crates/msa/tests/cases/projection_properties.rs)、[`test_msa_buffers.py`](../../src/tests/query/test_msa_buffers.py)。

## 读取原始记录

需要块的完整信息（`score`、`e` 行状态、`i` 行、负链坐标等）时，用 `IntervalReader` 读原始记录：

```python
from omicsweave.query import IntervalReader

with IntervalReader("multiz100way.bigmaf") as f:
    hits = f.fetch(("chr22", 10_600_000, 10_600_100), clip=False)
    blocks = [b.split(";") for b in hits.annotation.to_list()]
```

`clip=False` 保留每个完整块的参考区间，便于解释 `mafBlock` 中的坐标；`clip=True` 只改返回的 `start` / `end`，不会同步裁剪该字符串中的比对行。

## 查询耗时与块布局

以下是原生 MSA 初版的历史测量，不是当前审核后源码的重新计时。方法、原始日志与 JSON 可定位到本地任务产物 [`msa-validation.md`](../../tmp/msa-tabix-hicpro/msa-validation.md)、[`msa-measurements.json`](../../tmp/msa-tabix-hicpro/msa-measurements.json) 与 [`msa-cache-measurements.json`](../../tmp/msa-tabix-hicpro/msa-cache-measurements.json)。这些产物不属于发行接口；本轮增加了行元数据校验并调整输出分配，当前性能应重新成对测量，不能从下表推断普遍加速。

历史测量是在 CPU 上对 UCSC hg38 chr22 100way MAF 开头约 120 MB 压缩数据转换得到的 bigMaf 实测：208,903 个块，参考范围 10,510,000–20,535,489，固定输出全部 100 个物种。文件位于 tmpfs；每种窗口选择 30 个随机基因组起点和 30 个随机比对块起点，预热一轮后测三轮均值。时间包括结果分配，单位为 ms/次；具体耗时随块大小、覆盖密度和存储变化。

| 窗口 | Python 投影基线，无缓存 | 初版 Rust 投影，无缓存 |
|---|---:|---:|
| 128 bp | 1.929 | 1.897 |
| 512 bp | 2.909 | 2.264 |
| 131,072 bp | 192.97 | 46.67 |

在同一组 128 bp 窗口上，使用可容纳其解码工作集的 2 GiB `BlockCache` 并预热后，基线为 0.742 ms，Rust 为 0.353 ms。缓存的是解码块，仍需为每次查询投影；缓存放不下工作集时，随机窗口的收益可能很小。

默认 bigBed 每块 512 条记录。小窗口往往只需要几条，却仍需读取和解压整个数据块。通过 `ow.convert.run(..., bbi=BbiOptions(items_per_slot=64))` 可以选择更细的块；同一份输入的试验如下，默认值保持 512：

| `items_per_slot` | bigMaf 字节数 | 相对默认大小 | 128 bp，无缓存 | 128 bp，128 MiB 缓存 |
|---:|---:|---:|---:|---:|
| 512（默认） | 123,631,555 | 1.000× | 1.897 | 1.773 |
| 64 | 128,301,114 | 1.038× | 0.548 | 0.355 |
| 16 | 144,243,094 | 1.167× | 0.428 | 0.316 |
| 1 | 438,522,930 | 3.547× | 0.520 | 0.363 |

更细的块减少小窗口解压量，同时降低压缩率、增加索引开销；单条记录一个块在此样本中没有比 16 条一个块更快。以上大小仅包含 bigMaf 文件，不含 sidecar。

## 相关格式

- [MAF](maf.md)：bigMaf 的来源格式。
- [bigBed](bigbed.md)：bigMaf 使用的容器。
- [bigWig](bigwig.md)：phyloP、phastCons 等由比对计算的保守性分数。

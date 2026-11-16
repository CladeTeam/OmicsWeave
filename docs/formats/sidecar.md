# Sidecar：可选的加载元信息

[格式总览](README.md) · [查询](../query.md) · [转换](../convert.md) · [bigWig](bigwig.md) · [bigBed](bigbed.md) · [.hic](hic.md)

`.sidecar` 是 OmicsWeave 自己定义的**伴随文件**，放在 bigWig、bigBed 或 `.hic` 文件旁边，记录源文件中每个数据块的位置、压缩大小、**解压后的大小**和记录条数。这些信息源格式本身并不保存，而自行调度 I/O 的数据加载器在读取任何数据块之前就需要它们，例如预先分配解压缓冲区、估算一个批次要读多少字节，或决定哪些块值得缓存。

sidecar 是**可选的**。普通的 `TrackReader`、`IntervalReader`、`ContactReader` 查询完全不需要它，也不会去找它；它也不能替代 FASTA、BAM、BCF、CRAM 所需的索引。不确定是否需要时，就不必生成。

## 速览

| 项目 | 说明 |
|---|---|
| 用途 | 描述 bigWig / bigBed / `.hic` 的数据块布局，供数据加载器提前规划 |
| 文件名 | `<源文件>.sidecar`，如 `signal.bw.sidecar` |
| 存储 | 小型 little-endian 二进制，固定宽度记录 |
| 适用的源格式 | bigWig、bigBed、`.hic` |
| 是否必需 | 否；普通区域查询不使用 |
| 与源文件的绑定 | 记录源文件长度和内容指纹，读取时校验 |
| 生成方式 | 转换时 `--sidecar`、writer 的 `sidecar=True`、quant 的 `--sidecar`，或对已有文件 `build_sidecar()` |

## 为什么需要它

BBI（bigWig / bigBed）和 `.hic` 都把数据切成 zlib 压缩块，并用索引告诉读取方“某个区域在哪几个块、每块从哪个偏移开始、占多少压缩字节”。但它们都**不记录块解压后有多大、里面有多少条记录**：

- BBI 的 R 树索引只给出块的文件偏移和压缩长度。header 中的 `uncompressBufSize` 只是全文件的最大值。
- `.hic` 的块索引只记录块的位置和压缩字节数；zlib 流本身不带解压长度，而且每个块的内部编码（按行列表或稠密矩阵、短字段或宽字段）是逐块选择的，不解压就无法知道它解出来多大。

对于一次查询一个区域的普通读取，这不是问题：读出块、解压、解码即可。但对于需要在大量窗口间共享解压结果、按字节预算管理缓存、或在后台线程中预取的加载器，提前知道每块的解压大小和条目数，可以避免先解压一遍再分配。sidecar 就是把这些信息**测量一次**并保存下来。

## 文件结构

sidecar 由固定的容器头和按源格式定义的正文组成，所有字段都是 little-endian，8 字节字段保持 8 字节对齐，便于直接按表读取。

### 容器头（48 字节）

| 偏移 | 类型 | 字段 |
|---|---|---|
| 0 | `[u8; 8]` | magic：`OWSIDECR` |
| 8 | `u16` | 容器版本 |
| 10 | `u16` | 正文类型：BBI 或 hic |
| 12 | `u32` | 正文版本 |
| 16 | `u64` | 源文件长度（字节） |
| 24 | `u64` | 源文件内容指纹 |
| 32 | `u64` | 正文长度 |
| 40 | `u64` | 正文校验和 |
| 48 | | 正文 |

### bigWig / bigBed 正文

按层组织：先是原始数据层，再依次是各 zoom 层。每层记录其 reduction（原始层为 0）和块数，随后是每个块的一条记录：

| 字段 | 含义 |
|---|---|
| 起始染色体编号、起点、结束染色体编号、右端 | 索引记录的块跨度 |
| 文件偏移 | 块在源文件中的起始位置 |
| 文件中的字节数 | 压缩后大小 |
| 解压后的字节数 | 未压缩文件则等于存储大小 |
| 条目数 | bigWig 信号条目、bigBed 记录或 zoom 汇总记录的数量 |

正文还保存 header 的 `uncompressBufSize`（未压缩文件为 0）。块按索引叶子的顺序排列，块在其层中的位置即其序号。

### `.hic` 正文

记录 `.hic` 版本（8 或 9）和每个接触矩阵，矩阵按（染色体 1、染色体 2、单位、分辨率）排序：

| 矩阵字段 | 含义 |
|---|---|
| 染色体编号 1、2 | header 中的顺序，`chr1 <= chr2` |
| 单位、分辨率 | BP 或 FRAG，bin 大小 |
| 块网格 | 每块每边的 bin 数、每行的块数 |
| 块范围 | 该矩阵的块在块表中的起点与数量 |

| 块字段 | 含义 |
|---|---|
| 块编号 | `.hic` 中的块编号 |
| 压缩字节数、文件位置 | 块在源文件中的位置与大小 |
| 解压后的字节数 | 块解压后的大小 |
| 记录数 | 块解码后的接触记录数 |

## 与源文件的绑定

sidecar 描述的是某一个具体文件的字节布局，源文件一旦改写，旧 sidecar 就失效。为此容器头记录了源文件的长度和**内容指纹**：指纹由文件开头和末尾各 64 KiB、以及中间均匀分布的 14 个 4 KiB 页计算得出。

- 指纹只取决于内容，不取决于路径或修改时间：源文件与 sidecar 一起复制或移动后仍然有效。
- 读取 sidecar 时会重新计算源文件的长度和指纹，不一致时抛出 `StaleSidecarError`。
- 指纹是抽样的：长度不变、且所有抽样字节都未改变的改写无法被发现。替换源文件后应主动重新生成 sidecar，不要依赖校验兜底。

正文另有校验和，用于发现 sidecar 本身的损坏。sidecar 是派生产物，随时可以删除并从源文件重建。

## OmicsWeave 中的支持

| 功能 | 接口 |
|---|---|
| 文件对象 | `Sidecar`；`Sidecar.beside(path)` 构造 `<path>.sidecar` 的路径对象，不生成文件 |
| 格式识别 | `formats.detect` 按 magic 识别为 `Sidecar` |
| 为已有文件生成 | `BigWig(...).build_sidecar()`、`BigBed(...).build_sidecar()`、`Hic(...).build_sidecar()`，会读取整个源文件 |
| 写入时生成 | `BigWigWriter`、`BigBedWriter`、`HicWriter`、`HicBuilder` 的 `sidecar=True`，边写边记录，无需再读一遍 |
| 转换时生成 | `convert.run(..., sidecar=True)` / `ow convert --sidecar`，适用于 bigBed、bigWig、hic 目标 |
| quant 输出 | `quant.Coverage(sidecar=True)`、`quant.Contacts(sidecar=True)`，CLI `ow quant --sidecar` |
| 读取 | `formats.read_sidecar(源文件)` → `BbiSidecar` 或 `HicSidecar` |

所有生成方式默认都**不**写 sidecar。

## 生成

转换时一并生成：

```sh
ow convert signal.bedGraph signal.bw --sidecar --chrom-sizes genome.chrom.sizes
ow convert contacts.mcool contacts.hic --sidecar
```

为已有文件单独生成：

```python
from omicsweave.formats import BigWig, Hic

sidecar = BigWig("signal.bw").build_sidecar()
print(sidecar.path)  # signal.bw.sidecar

Hic("contacts.hic").build_sidecar()
```

`build_sidecar()` 需要读取并解压整个源文件来测量每个块，大文件耗时与完整读取一遍相当。如果文件是自己写出的，在 writer 上设置 `sidecar=True` 更省事：

```python
from omicsweave.formats import BigWigWriter

with BigWigWriter("signal.bw", {"chr1": 1000}, sidecar=True) as writer:
    writer.add("chr1", [100, 200], [150, 300], [2.5, 8.0])
# 离开 with 时同时写出 signal.bw 和 signal.bw.sidecar
```

## 读取

```python
from omicsweave.formats import read_sidecar

metadata = read_sidecar("signal.bw")  # 传源文件路径，而不是 .sidecar 路径
print(metadata.uncompress_buf_size)
for level in metadata.levels:          # 原始层在前，然后是各 zoom 层
    blocks = level.blocks
    print(level.reduction, len(level))
    print(blocks["offset"], blocks["size"], blocks["inflated"], blocks["items"])
```

`BbiSidecar.levels` 中的每个 `BbiLevelBlocks` 提供 `reduction` 和按列存放的 `blocks`，列名为 `start_chrom`、`start`、`end_chrom`、`end`、`offset`、`size`、`inflated`、`items`，都是 `uint64` 数组。

```python
metadata = read_sidecar("contacts.hic")
print(metadata.version)
for m in metadata.maps:
    print(m.chrom1, m.chrom2, m.unit, m.resolution, m.block_bin_count, m.block_column_count)
    print(m.blocks["number"], m.blocks["position"], m.blocks["size"], m.blocks["inflated"], m.blocks["records"])
```

`HicSidecar.maps` 中的每个 `HicMapBlocks` 描述一个矩阵；染色体用 header 中的编号表示，可与 `ContactReader(...).chroms` 的顺序对应。

例如，估算一个 bigWig 原始层全部解压后需要多少内存：

```python
level0 = read_sidecar("signal.bw").levels[0]
print(int(level0.blocks["inflated"].sum()), "bytes")
```

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 查询需要先生成 sidecar 吗 | 不需要；普通查询从不读取 sidecar |
| `StaleSidecarError` | sidecar 是为另一个版本的源文件生成的；重新 `build_sidecar()` |
| `read_sidecar` 报文件不存在 | 源文件旁没有 `<源文件>.sidecar`；先生成 |
| 传入 `.sidecar` 路径报错 | `read_sidecar` 接受的是源文件路径 |
| 能否给 BAM / FASTA 生成 | 不能；只支持 bigWig、bigBed、`.hic`，其他格式使用各自的索引 |
| 复制文件后是否还有效 | 与源文件一起复制或移动仍有效；指纹只看内容 |

## 相关格式

- [bigWig](bigwig.md)、[bigBed](bigbed.md)：BBI 容器，sidecar 记录其原始层与 zoom 层的块。
- [.hic](hic.md)：sidecar 记录每个接触矩阵的块。
- [BAI](bai.md)、[CSI](csi.md)、[FAI](fai.md)、[GZI](gzi.md)、[CRAI](crai.md)：其他格式的区域索引，与 sidecar 用途不同。

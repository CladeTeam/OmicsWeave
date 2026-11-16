# BAM

[格式总览](README.md) · [查询](../query.md) · [转换](../convert.md) · [定量](../quant.md) · [SAM](sam.md) · [BAI](bai.md) · [CSI](csi.md)

BAM 是 [SAM](sam.md) 比对记录的**压缩二进制形式**：字段与 SAM 一一对应，但以紧凑的二进制编码保存，并用 BGZF 分块压缩。按坐标排序并建立索引后，可以直接跳到某个区域，只解压包含该区域读段的少数几个块。它是测序比对数据最通用的交换和存储格式。

在 OmicsWeave 中，BAM 有两种用法：`AlignmentReader` 通过 BAI / CSI 索引**按区域读取**比对记录，返回列式 NumPy 数组；[quant](../quant.md) 则**顺序扫描**坐标排序 BAM，一次扫描同时生成覆盖度、计数矩阵等多种产物。SAM、CRAM 都先转换为 BAM 再使用。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 读段比对记录，字段同 SAM |
| 常见后缀 | `.bam` |
| 存储 | 二进制记录，BGZF 分块压缩 |
| 坐标 | 二进制中为 0-based；OmicsWeave 结果为 0-based 半开 |
| 随机访问 | 坐标排序 + 索引后支持 |
| 外部索引 | [BAI](bai.md)（`.bam.bai`）或 [CSI](csi.md)（`.bam.csi`） |
| 典型来源 | `samtools sort`、比对流程（STAR、Cell Ranger、bwa + samtools）、`ow convert` |
| OmicsWeave | 区域查询、建索引、SAM / CRAM → BAM、BAM 坐标排序、quant 输入 |

## 格式由来

BAM 与 SAM 同时由 samtools 引入（Li et al., 2009），规范同样维护在 hts-specs 中。它的核心设计是 **BGZF**（Blocked GNU Zip Format）：把数据切成不超过 64 KiB 的块，每块是一个独立的 gzip member。这样文件仍是合法的 gzip（`zcat` 可以解压），又能从任意块开始解压，配合索引实现随机访问。BGZF 后来也被 [BCF](bcf.md)、bgzip [FASTA](fasta.md)、tabix 索引的文本文件采用。

## 文件结构

### BGZF 分块压缩

```text
┌─────────────┬─────────────┬─────────────┬─────┬──────────────┐
│ BGZF 块 0   │ BGZF 块 1   │ BGZF 块 2   │ ... │ EOF 空块     │
│ ≤64 KiB 解压 │             │             │     │ 28 字节      │
└─────────────┴─────────────┴─────────────┴─────┴──────────────┘
```

- 每个块是带 `BC` 扩展字段的 gzip member，扩展字段记录本块压缩后的大小，因此不解压也能跳到下一个块。
- 文件以一个固定的 28 字节空块结尾，用来判断文件是否完整。
- 块彼此独立，所以可以多线程并行压缩和解压。

**虚拟偏移**（virtual offset）是 BAM 中定位一条记录的方式：

```text
virtual_offset = (压缩块在文件中的字节偏移 << 16) | 记录在解压块内的偏移
```

高 48 位找到块，低 16 位找到块内位置。BAI / CSI 中保存的都是虚拟偏移。

### 解压后的数据流

```text
magic "BAM\1"
l_text + header 文本（与 SAM header 相同）
n_ref + 每条参考序列的名称和长度
记录 1
记录 2
...
```

### 一条记录的二进制布局

| 字段 | 类型 | 对应 SAM 列 / 说明 |
|---|---|---|
| `block_size` | `int32` | 本记录剩余字节数 |
| `refID` | `int32` | RNAME 的参考编号（tid），`-1` 表示 `*` |
| `pos` | `int32` | POS − 1，0-based |
| `l_read_name` | `uint8` | QNAME 长度（含结尾 NUL） |
| `mapq` | `uint8` | MAPQ |
| `bin` | `uint16` | 该记录所属的 BAI bin |
| `n_cigar_op` | `uint16` | CIGAR 操作数 |
| `flag` | `uint16` | FLAG |
| `l_seq` | `int32` | 序列长度 |
| `next_refID` | `int32` | RNEXT 的参考编号 |
| `next_pos` | `int32` | PNEXT − 1 |
| `tlen` | `int32` | TLEN |
| `read_name` | `char[]` | QNAME，NUL 结尾 |
| `cigar` | `uint32[]` | 每个操作编码为 `长度 << 4 \| 操作码`，操作码按 `MIDNSHP=X` 排列 |
| `seq` | `uint8[]` | 每个碱基 4 bit，两个碱基一个字节 |
| `qual` | `uint8[]` | 原始 Phred 值（不加 33）；缺失时全部为 `0xFF` |
| aux | 字节 | 可选标签：两字节标签名 + 一字节类型 + 值 |

`n_cigar_op` 只有 16 位。CIGAR 超过 65535 个操作时（长读长常见），真实 CIGAR 保存在 `CG` 标签中，记录本身放一个占位 CIGAR；OmicsWeave 读取时会自动用 `CG` 替换占位符。

## 关键概念

### 排序与索引

只有**按坐标排序**的 BAM 才能建索引：记录先按 tid（`@SQ` 顺序）排列，tid 为 `-1` 的未定位读段放在最后，同一参考内按位置递增。索引记录每个基因组 bin 的记录从哪个虚拟偏移开始，查询时只需读取相关块。

比对软件直接输出的 BAM 通常未排序或按读段名排序，需要先排序：

```sh
ow convert aligned.bam sorted.bam --sort --index -t 8 -m 4G
```

### 重叠查询的语义

区域查询返回**与区域重叠**的记录：记录的参考跨度 `[start, end)` 与查询区域有交集。跨度由 POS 和 CIGAR 中消耗参考的操作（`M`、`D`、`N`、`=`、`X`）决定，因此跨越内含子（`N`）的 RNA-seq 读段，即使查询区域落在内含子中也会返回。没有 CIGAR 的记录跨度按 1 bp 计算。

放在 mate 位置上的未比对读段（FLAG 含 `0x4` 但有 RNAME / POS）同样会被返回；RNAME 为 `*` 的读段无法按区域取到。

### 过滤是使用者的事

BAM 中通常同时存在主比对、次要比对、补充比对、重复读段和低质量比对。区域查询**原样返回**所有记录，与 `samtools view region` 一致；需要什么样的读段，由使用者根据 FLAG、MAPQ 等列筛选，或交给 quant 的过滤规则。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 格式识别 | 解压首块后按 `BAM\1` magic 识别为 `Bam` |
| 区域查询 | `AlignmentReader.fetch`，返回记录与 htslib 区域迭代器（`sam_itr_queryi`）一致，按文件顺序 |
| 列式结果 | 固定列为 NumPy 数组；名称、序列、质量、CIGAR、标签按需解码 |
| 建索引 | `Bam.build_index()`：BAI 或 CSI，与 `samtools index` 的 bins / chunks 相同 |
| 索引查找 | 默认依次找 `<bam>.csi`、`<去扩展名>.csi`、`<bam>.bai`、`<去扩展名>.bai`；可用 `index=` 指定 |
| header | `chroms`、`references`、`header_text` |
| 缓存与并发 | 可接 `BlockCache` 共享解压后的 BGZF 块；解码释放 GIL；可 `clone()`、pickle |
| 转换 | SAM → BAM、CRAM → BAM、BAM → 坐标排序 BAM，均可同时建索引 |
| quant | 坐标排序 BAM 作为输入，顺序扫描，不需要索引 |

限制：

- 不写 SAM，不写 CRAM；没有 BAM → SAM / CRAM 转换。
- BAM → BAM 只用于排序，必须指定 `sort=True` / `--sort`；不提供纯复制或按读段名排序。
- 打开 reader 不会自动建索引；缺少索引时报 `IoError`，提示调用 `build_index`。

## 准备：排序与索引

```sh
# 未排序 BAM：排序并建索引
ow convert unsorted.bam sorted.bam --sort --index

# 指定 CSI 及其最小 bin
ow convert unsorted.bam sorted.bam --sort --index csi --min-shift 14
```

已有坐标排序 BAM，只缺索引：

```python
from omicsweave.formats import Bam

index = Bam("sorted.bam").build_index()       # kind="auto"
print(index)                                  # Bai('sorted.bam.bai')，参考过长时为 Csi
Bam("sorted.bam").build_index("csi", min_shift=14, threads=4)
```

`kind="auto"` 在每条参考序列都不超过 `2**29` bp 时写 BAI，否则写 CSI。`threads` 默认使用进程可用的全部 CPU 核心来解压文件。输入未排序时报错。

排序行为与 `samtools sort` 相同：按参考编号（未定位读段最后）、位置、正链先于负链排列，相同键保持输入顺序；header 改为 `SO:coordinate`，不添加 `@PG`。超出 `memory` 预算时把已排序的批次写入临时文件再归并，输出与线程数、是否溢写无关。临时文件默认放在输出旁边，可用 `--tmp-dir` / `tmp_dir` 指定。

## 区域查询

```python
from omicsweave.query import AlignmentReader

with AlignmentReader("sorted.bam") as reads:
    print(reads.chroms)        # {'chr1': 248956422, ...}，按参考编号顺序
    print(reads.references)    # ('chr1', 'chr2', ...)，mate_tid 的下标
    hits = reads.fetch(("chr1", 100, 200), fields=("name", "seq", "cigar"))
    print(len(hits), hits.start, hits.end, hits.flag, hits.mapq)
    if len(hits):
        print(hits.names[0], hits.seqs[0], hits.cigars.string(0))
```

### 固定列

默认只解码固定列，速度最快：

| 字段 | dtype | 含义 |
|---|---|---|
| `start` / `end` | `int64` | 比对起点与参考跨度右端，0-based 半开 |
| `flag` | `uint16` | SAM FLAG，可与 `types.SamFlag` 做位运算 |
| `mapq` | `uint8` | 比对质量 |
| `mate_tid` | `int32` | mate 的参考编号，对应 `reads.references`；`-1` 表示没有 |
| `mate_start` | `int32` | mate 的 0-based 起点；`-1` 表示没有 |
| `tlen` | `int32` | SAM TLEN |

### 可变长字段

通过 `fields` 选择，可以是名称元组、`types.AlignmentField` 标志组合，或 `"all"`：

| 名称 | 结果属性 | 类型与内容 |
|---|---|---|
| `name` | `names` | `StringVector`，读段名 |
| `seq` | `seqs` | `StringVector`，ASCII 碱基 |
| `qual` | `quals` | `StringVector`，原始 Phred 字节（不是 FASTQ 的 +33 编码）；缺失为 `0xff` |
| `cigar` | `cigars` | `Cigars`：BAM 编码的 `uint32` 操作与 `offsets` |
| `aux` | `aux` | `StringVector`，BAM 二进制编码的可选标签 |

未请求的字段为 `None`。

```python
import numpy as np
from omicsweave.types import AlignmentField, SamFlag

with AlignmentReader("sorted.bam") as reads:
    hits = reads.fetch(("chr1", 0, 100_000), AlignmentField.QUAL | AlignmentField.CIGAR)

    # 主比对、非重复、MAPQ ≥ 30
    bad = SamFlag.UNMAPPED | SamFlag.SECONDARY | SamFlag.SUPPLEMENTARY | SamFlag.DUPLICATE
    keep = ((hits.flag & bad) == 0) & (hits.mapq >= 30)

    phred = np.frombuffer(hits.quals.bytes(0), dtype=np.uint8)  # 第 0 条的质量
    ops = hits.cigars[0]                                        # uint32 数组
    lengths, codes = ops >> 4, ops & 0xF                        # 操作码按 "MIDNSHP=X"
    print(hits.cigars.string(0))                                # 例如 "3S5M2D4M"
```

aux 是 BAM 的二进制标签编码，例如 `NM:i:0` 存为 `b"NMC\x00"`（标签名、类型 `C` 表示 uint8、值 0），需要按 SAM 规范解析。

### 链方向与顺序

记录按文件顺序（即坐标顺序）返回。负链 key 不会改变记录方向，也不会反向互补 SEQ；SEQ 始终是参考正链方向，读段本身的方向看 FLAG 的 `0x10`。

### 自定义索引位置

```python
with AlignmentReader("sorted.bam", index="index/sorted.csi") as reads:
    print(reads.index)  # Csi('index/sorted.csi')
```

## 定量

quant 接受一个或多个坐标排序 BAM，顺序扫描，不要求 BAI / CSI。覆盖度、RNA / ATAC / ChIP 计数、剪接、甲基化、Hi-C 和等位计数都在一次扫描中完成：

```sh
ow quant sorted.bam -o results --coverage --gtf genes.gtf --genes bulk
```

各产物有自己的默认过滤（例如覆盖度默认排除 unmapped、secondary、QC-failed、duplicate），见 [quant](../quant.md)。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 提示缺少索引 | 调用 `Bam(path).build_index()`；reader 不会自动建索引 |
| 建索引报未排序 | 先 `ow convert in.bam sorted.bam --sort` |
| 找不到染色体 | 查看 `reads.chroms`，使用 header 中的名称（`chr1` 与 `1` 不互换） |
| 结果中有重复和次要比对 | 查询不过滤；用 `flag` / `mapq` 列筛选 |
| 质量值看起来比 FASTQ 小 33 | `quals` 是原始 Phred 值，不是 ASCII 编码 |
| 负链读段序列不是读段原始方向 | BAM 中 SEQ 总是参考正链方向；需要原始方向时按 `0x10` 自行反向互补 |
| 参考序列超过 512 Mb 无法写 BAI | 使用 CSI：`--index csi` |

## 相关格式

- [SAM](sam.md)：BAM 的文本形式，字段含义的完整说明。
- [CRAM](cram.md)：参考压缩的比对格式，可转为 BAM。
- [BAI](bai.md)、[CSI](csi.md)：BAM 的区域索引。
- [bigWig](bigwig.md)：quant 覆盖度的输出格式。

# 基因注释

[文档首页](README.md) · [查询](queries.md) · [定量分析](quantification.md)

## 从注释读取转录本

```python
from omicsweave import Region
from omicsweave.annotation import Annotation
from omicsweave.query import SequenceReader

ann = Annotation.read("genes.gtf.gz", attributes="all", threads=4)
print(len(ann.genes), len(ann.transcripts), ann.report.problems)
coding = ann.subset(biotypes=["protein_coding"], tags=["basic"])
keys = ann.loci(coding.transcripts, "exons")
with SequenceReader("genome.fa.gz") as genome:
    if keys:
        sequence = genome.fetch(keys[0]).text
rows = ann.overlaps(Region("chr1", 1000, 2000), level="exon")
starts = ann.exons.start[rows]
```

`Annotation.read` 完整载入源文件，并建立可复用的列式模型；它不按查询窗口重新读取 GTF。返回的坐标统一为 0-based、半开区间。`overlaps` 返回所选表的行号，不是特征对象；`genes.id` 等字符串列是 `StringVector`，可用 `.to_list()`。

| 需求 | 接口 |
|---|---|
| id / 名称查找 | `gene(key)`、`transcript(key)` 返回第一个匹配；用 `gene_rows` / `transcript_rows` 获取全部匹配 |
| 区间查找 | `overlaps`、`overlaps_batch`、`nearest` |
| 筛选 | `subset(biotypes=..., tags=...)`，返回行号集合 |
| 转录本序列 key | `locus(transcript, part)`、`loci(rows, part)`；`part` 为 `exons` / `cds` / `utr5` / `utr3` |
| 派生区域 | `introns`、`exonic_regions`、`promoters`、`gene_bodies`、`tss`、`tes`、`tes_windows` |
| 剪接与外显子 bin | `junctions`、`exon_bins` |
| 染色体名对齐 | `match_chroms`；可传别名表 |
| 导出 | `write_gtf`、`write_bed12` |

只做定量或区间查询时可用 `attributes=None` 降低内存；id、name、biotype、tag 仍保留。需要原样写回 GTF 时使用 `attributes="all"` 和默认的原始布局。`strict=False` 默认收集结构问题到 `report`；语法错误始终报错。`strict=True` 遇到结构问题也会报错。

`Annotation` pickle 后在目标进程按路径重新读取；其内存表不随 pickle 共享。

## 解析与数据模型

`omicsweave-gtf` 把基因注释（GTF 2.2、GFF3）读成一个通用的、可复用的内存模型，供定量（BAM 读段按基因 / 外显子 / 内含子 / 剪接位点计数）、单细胞 ATAC 的 gene activity、按注释采样等功能共用。它不是训练数据加载的一部分：注释整体常驻，按区间和按 id 查询。

### 读取

- 支持 GENCODE、Ensembl、RefSeq、UCSC 导出的 GTF 和 GFF3；纯文本、gzip、BGZF 按内容识别。方言按 `##gff-version 3` 或第一条记录的属性语法识别，也可显式指定。
- 文本按行切批，在 `omicsweave_core::exec` 的线程池上并行解析，每批得到使用批内字符串编号的列；合并在单独的线程上按批次顺序进行，把批内字符串依次并入全局字符串池。所以字符串编号、行顺序与线程数和批大小无关，测试检查这一点。普通 gzip 在单独线程上解压，BGZF 在池上并行解压。
- 规范化：坐标统一为 0-based、左闭右开的 `u32`（只在读写时换算）；链为 `Strand`（`+`、`-`、`.`、`?`）；phase 为 0–2 或无；score 按原文保存（`.` 为空），因为 UCSC 写成 `0.000000`，转成浮点就不能原样写回。
- 层级是 gene → transcript → feature。GTF 按 (seqid, `gene_id`) 和 (seqid, `transcript_id`) 分组，同一个 id 出现在两条染色体上（PAR、UCSC 的多位置映射）就是两个基因；RefSeq 基因行的 `transcript_id ""` 视为没有。GFF3 按 `ID` / `Parent` 建树，转录本靠结构而不是类型名识别：有外显子级子记录（exon、CDS、UTR、密码子等）的记录是转录本（`mRNA`、`lnc_RNA`、`V_gene_segment` 都可以），转录本的父记录是基因；同一 `ID` 的多行外显子级记录（RefSeq 的 CDS）是一个不连续特征的多行；有多个父转录本的特征在每个转录本下各有一行。
- 缺失的 gene / transcript 记录由子记录合成：范围取子记录的最小 / 最大坐标，链取第一个子记录的链（UCSC 导出没有 gene / transcript 行；RefSeq 线粒体基因的 CDS 直接挂在 gene 下，合成一个转录本）。
- 名称、biotype、tag 从各来源的常用键读取：`gene_name` / `gene` / `Name`，`gene_type` / `gene_biotype` / `biotype`，`transcript_*` 同理；`tag` 可以重复（GTF）或逗号分隔（GFF3）。常见 tag（`basic`、`Ensembl_canonical`、`MANE_Select`、`CCDS`、`appris_principal_*`、`RefSeq Select` 等）同时存为位图，大小写和空格 / 下划线的差别不影响匹配，所有 tag 值另以字符串列表保存。这些字段与属性保留策略无关，总会读取。
- 染色体别名：读取时可按表改名；`match_seqids` 把另一个文件的名字（BAM header 的参考名）对应到注释的 seqid，依次试原名、别名表（双向）、内置的 UCSC ↔ Ensembl 写法（`chr1` ↔ `1`、`chrM` ↔ `MT`）。
- 校验默认只收集报告：孤儿记录、未知 `Parent`、重复 id、链冲突、seqid 冲突、子记录超出父记录、不合层级的嵌套、多个父转录本。`strict` 时第一个问题即报 `Error::Format`。输入不要求排序。语法错误（列数、坐标、链、phase、引号）总是报错。

### 内存模型

不为每条记录分配堆对象：
- `StringPool`：一块字节 arena 加偏移，内容查找用开放寻址的 id 表。seqid、source、类型、id、名称、biotype、属性键和值都只存一次，用 `u32` 引用。
- `Genes`、`Transcripts`、`Features`、`Others` 是列式结构。一个基因的转录本、一个转录本的特征都是连续的行区间（CSR 偏移），按位置排序；基因按 (seqid, start, end) 排序，所以基因行号即基因组顺序。不属于层级的记录（GFF3 的 `region`、`chromosome`，有问题的记录）进 `Others`，保留文件顺序。
- 属性（第 9 列）按文件记录编号存成一张平表：每个不同的 (key, value) 对只存一次，记录只存对的编号。GENCODE 外显子行重复转录本的十几个属性，按对去重后，1.6M 行、全部属性的模型比直接存 (key, value) 小约 30%。保留策略 `All`（无损、原顺序，写回原 GTF 需要）、`Only(keys)`、`None`。GTF 值是否带引号记在键编号的最高位。
- 按 id、按名称的查找表首次使用时构建；索引、派生视图同样按需构建并缓存。注释构建完成后不可变，可在线程间共享。

### 区间索引

`IntervalIndex` 按 seqid 分组，组内按 (start, end, id) 排序存四个平行数组（start、end、id、子树最大 end），即 cgranges 的隐式增强区间树：树就是数组本身，构建是一次排序加一次自底向上的扫描，查询只访问连续内存。基因、转录本、外显子各有一个索引，id 是对应表的行号。查询：overlap（按 start 顺序给出）、contained_in、containing、nearest（重叠为 0，否则取两侧最近的全部并列项）、批量查询（CSR 输出）。结果追加到调用方的缓冲区。

### SweepCursor

坐标排序的 BAM 逐条扫描时，读段起点不减。游标维护「已进入」（起点小于某次查询终点）且「未过期」（终点大于下限 floor）的活动集合：每个项目只进入、过期各一次，推进是摊还 O(1)；一次查询的代价是活动集合的大小，不分配内存。剪接读段的多个比对块彼此不单调，所以分成两步：每条读段 `advance(read_start)` 一次，再对每个块调用 `overlaps`；起点不低于 floor 的任何查询都精确。低于 floor 的查询从树重新定位（O(log n + k)），结果仍然正确。

### 派生视图

按需计算，无参数的视图缓存在注释上：
- 每个转录本的外显子块（相交或相接的外显子合并；没有 exon 行时用 CDS / UTR / 密码子行），内含子，基因的合并外显子区，注释剪接位点集合（每个不同的 (seqid, start, end, strand) 一次，附 donor / acceptor、所属基因和转录本，可按坐标精确查找），DEXSeq 式不相交外显子 bin（按基因，或按同链外显子重叠的基因聚合），转录本长度、CDS 长度。
- 计数特征窗口：TSS / TES，promoter（`GenomicRanges::promoters` 的定义）和 TES 窗口（上下游可配置、按链方向），基因体加上游延伸（ArchR / Signac 的 gene activity）；在 0 处截断，给出染色体长度时在末端截断。没有链信息的记录按 `+` 处理。
- 区间视图都是 `IntervalSet`：(seqid, start, end, strand, owner, name) 列加自己的索引和游标，owner 是所属基因或转录本的行号，所以计数层对它们和对基因、外显子用完全相同的方式索引和扫描。peak、固定 bin 等其它来源的特征也用这个类型（`IntervalSetBuilder`、`IntervalSet::from_columns`）。
- `Derive` 在 `Selection` 上、按染色体长度截断地构建这些视图。`select` 按 biotype、tag、基因 id / 名称、染色体给出基因和转录本的行号列表，不复制数据。

### 输出

GTF 有两种布局：`Original` 按文件顺序写回全部记录和 `#` 行，对全部属性读入的 GTF 逐字节一致（合成的记录不写）；`Model` 按模型写 gene、transcript、feature 行，先写模型字段（`gene_id`、`transcript_id`、名称、biotype、tag），再写其余保留的属性，GFF3 也由此转成 GTF。BED12 按转录本写，block 是外显子块，thick 部分覆盖 CDS 和密码子，非编码转录本按 UCSC 惯例 `thickStart = thickEnd = chromEnd`。

### Python

`omicsweave.annotation.Annotation.read(path, attributes="all" | [keys] | None, aliases=None, strict=False)`：`genes`、`transcripts`、`features`、`exons` 是列式 dataclass（numpy 列，字符串为 `StringVector`），`gene(id)`、`overlaps(region, level)`、`nearest`、`subset(...)`、`introns()`、`junctions()`、`exon_bins()`、`promoters()`、`gene_bodies()`、`match_chroms()`、`report`。原生层只收发基本类型和 numpy 数组；注释不可变，不需要锁，调用时释放 GIL。pickle 时按路径和选项重新读取。

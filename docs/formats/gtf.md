# GTF

[格式总览](README.md) · [查询](../query.md) · [定量](../quant.md) · [GFF3](gff3.md) · [bigBed](bigbed.md)

GTF（Gene Transfer Format）是描述**基因结构**的九列文本格式：每行一个特征（基因、转录本、外显子、CDS、UTR、起止密码子……），用 `gene_id`、`transcript_id` 等属性把它们串成“基因 → 转录本 → 外显子”的层级。GENCODE、Ensembl、RefSeq、UCSC 都以 GTF 发布注释，STAR、featureCounts、Cell Ranger 等工具也以它为输入。

在 OmicsWeave 中，GTF 由 `Annotation.read` 一次载入内存，成为列式、带索引的注释模型。它不是按区域读取的数据文件，而是**生成查询 key 和定量特征的来源**：按转录本得到拼接外显子的 `Spliced` key 去取序列和信号，生成内含子、启动子、基因体等区间，或交给 `quant` 构建基因计数矩阵。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 基因、转录本及其子特征的位置和属性 |
| 常见后缀 | `.gtf`、`.gtf.gz` |
| 存储 | tab 分隔文本，九列 |
| 坐标 | 文件中 1-based、两端闭合；载入后为 0-based 半开 |
| 随机访问 | 不需要；整份载入内存后建索引 |
| 外部索引 | 不需要 |
| 典型来源 | GENCODE、Ensembl、RefSeq、UCSC、Cell Ranger 参考包 |
| OmicsWeave | `Annotation.read` 载入；查找、区间推导、查询 key、quant 特征；导出 GTF / BED12 |

## 格式由来

GTF 源自 GFF（General Feature Format）第二版，由 Sanger / WashU 的基因预测社区约定，用于交换基因模型。当前通行的是 GTF 2.2：九列结构沿用 GFF2，第九列属性必须包含 `gene_id` 和 `transcript_id`，以此表达层级。后来的 [GFF3](gff3.md) 改用 `ID` / `Parent` 显式表达层级，两者描述的信息大体相同。

GTF 2.2 规范本身只要求 `exon`、`CDS`、`start_codon`、`stop_codon` 等行，基因和转录本的范围由子特征推出；GENCODE / Ensembl 的 GTF 额外写出了 `gene` 和 `transcript` 行。

## 文件结构

### 九列

```text
chr17	HAVANA	gene	7661779	7687538	.	-	.	gene_id "ENSG00000141510.18"; gene_type "protein_coding"; gene_name "TP53";
chr17	HAVANA	transcript	7661779	7687538	.	-	.	gene_id "ENSG00000141510.18"; transcript_id "ENST00000269305.9"; gene_name "TP53"; transcript_type "protein_coding"; tag "basic"; tag "MANE_Select";
chr17	HAVANA	exon	7687377	7687538	.	-	.	gene_id "ENSG00000141510.18"; transcript_id "ENST00000269305.9"; exon_number 1;
chr17	HAVANA	CDS	7676521	7676594	.	-	0	gene_id "ENSG00000141510.18"; transcript_id "ENST00000269305.9"; exon_number 4;
```

| 列 | 名称 | 含义 |
|---|---|---|
| 1 | seqname | 染色体或 contig 名，需与参考基因组一致 |
| 2 | source | 注释来源，如 `HAVANA`、`ENSEMBL`、`BestRefSeq` |
| 3 | feature | 特征类型：`gene`、`transcript`、`exon`、`CDS`、`UTR`、`start_codon`、`stop_codon`、`Selenocysteine` 等 |
| 4 | start | 起点，1-based，包含 |
| 5 | end | 终点，1-based，包含 |
| 6 | score | 分值，通常为 `.` |
| 7 | strand | `+`、`-`，或 `.` 表示无方向 |
| 8 | frame | CDS 的读码相位 `0` / `1` / `2`：从该行起点跳过几个碱基到达下一个完整密码子；非 CDS 为 `.` |
| 9 | attribute | `key "value";` 形式的属性，用分号分隔 |

### 属性语法

```text
gene_id "G1"; transcript_id "T1"; exon_number 2; tag "basic"; tag "CCDS"; note "a; b";
```

- 每个属性是 `key value;`，值通常加双引号，数字值可以不加引号。
- 同一个 key 可以出现多次，例如多个 `tag`。
- 引号内可以包含分号。
- 不同来源的 key 名称不同：基因类型在 GENCODE 中是 `gene_type`，在 Ensembl 中是 `gene_biotype`；基因名可能是 `gene_name` 或 `gene`。

### 层级

```text
gene  (gene_id)
└── transcript  (transcript_id)
    ├── exon
    ├── CDS        编码部分，GENCODE / Ensembl 中不含终止密码子
    ├── UTR / five_prime_UTR / three_prime_UTR
    ├── start_codon
    └── stop_codon
```

同一转录本的所有行共享 `transcript_id`，同一基因的转录本共享 `gene_id`。负链基因的外显子在文件中可以按基因组坐标或按转录方向排列，阅读时要以 strand 判断 5′ / 3′。

### 坐标换算

GTF 的 `start=1001, end=2000` 表示第 1001 到第 2000 个碱基（共 1000 个），换成 0-based 半开坐标就是 `[1000, 2000)`。OmicsWeave 载入时统一换算，之后所有 start / end 都是 0-based 半开。

## 关键概念

### 基因、转录本、外显子的“部分”

一个转录本可以切成几种拼接区间：

| part | 含义 |
|---|---|
| `exons` | 全部外显子拼接：成熟转录本 |
| `cds` | 编码区块，**包含终止密码子**，因此完整 ORF 的长度是 3 的倍数 |
| `utr5` / `utr3` | 编码区 5′ / 3′ 侧的外显子碱基 |

GTF 中 `CDS` 行通常不含终止密码子（单独写为 `stop_codon`），OmicsWeave 生成 `cds` key 时会把它并入。

### 推导区间

很多分析需要的区间不在文件中直接写出：内含子是相邻外显子之间的间隔，启动子是 TSS 上下游的窗口，基因体加上游扩展用于 scATAC 的 gene activity。OmicsWeave 从载入的模型按链方向计算这些区间，返回 `IntervalSet`。

### 转录本标签

GENCODE / Ensembl / RefSeq 用 `tag` 标出代表性转录本，例如 `basic`、`Ensembl_canonical`、`MANE_Select`、`CCDS`。按标签筛选可以得到“每个基因一个转录本”的集合。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 读取 | `Annotation.read`：GENCODE、Ensembl、RefSeq、UCSC 的 GTF；纯文本、gzip、BGZF，按内容判断 |
| 层级 | 按 (seqid, `gene_id`) 和 (seqid, `transcript_id`) 组装；没有自身记录的基因 / 转录本由子特征合成 |
| 并行解析 | 文本切成批次多线程解析，结果与线程数无关 |
| 属性 | 默认全部保留；可只保留指定 key 或不保留 |
| 查找 | 按基因 / 转录本 id 或名称；按区域重叠、最近邻 |
| 查询 key | `locus` / `loci` 生成 `Spliced`，直接交给序列、信号 reader |
| 推导区间 | 内含子、外显子并集、剪接位点、DEXSeq 式 exon bins、启动子、TSS / TES、基因体 |
| 筛选 | 按 biotype、转录本 biotype、标签、基因、染色体 |
| 染色体名 | 载入时改名；与 BAM header 匹配时内置 `chr1` ↔ `1`、`chrM` ↔ `MT` |
| 质量检查 | 孤立记录、重复 id、链冲突、超出父区间等写入 `report`，或 `strict=True` 报错 |
| 导出 | `write_gtf`（原始记录顺序或按模型）、`write_bed12` |
| quant | 基因外显子、基因跨度、启动子、gene activity、velocity、ASE 等特征 |

GTF 不通过 `formats.detect` 或 `query.open` 打开（内容会被当作普通文本表），也不能交给区域 reader；始终使用 `Annotation.read`。

## 载入

```python
from omicsweave.annotation import Annotation

annotation = Annotation.read("gencode.v44.annotation.gtf.gz")
print(annotation)              # Annotation('...', dialect='gtf', genes=..., transcripts=...)
print(annotation.chroms[:3])   # 染色体名，按首次出现的顺序
print(annotation.heap_bytes)   # 注释模型占用的内存
```

`Annotation.read` 的参数：

| 参数 | 默认 | 含义 |
|---|---|---|
| `attributes` | `"all"` | `"all"` 保留全部属性（`original` 布局需要）；属性名列表只保留这些；`None` 不保留。id、名称、biotype、tag 总会读取 |
| `aliases` | `None` | 载入时改染色体名，例如 `{"1": "chr1"}`；`read_aliases(path)` 可从两列文件读取 |
| `strict` | `False` | 遇到第一个结构问题即抛出 `FormatError` |
| `threads` | 可用核心数，最多 8 | 并行解析线程数 |
| `dialect` | 自动 | 强制按 `"gtf"` 或 `"gff3"` 解析 |
| `batch_bytes` | 4 MiB | 每个解析批次的文本量 |

只做定量或生成 key 时，`attributes=None` 能明显减少内存。`Annotation` 不可变，可以在线程间共享；pickle 时只保存路径和读取参数，在 DataLoader worker 中会重新读取文件。fork 之后应在子进程中重新载入（或通过 pickle 传入）。

### 检查读取结果

```python
report = annotation.report
print(report.records, report.genes, report.transcripts, report.features)
print(report.synthesized_genes, report.synthesized_transcripts)
print(report.is_clean)
for problem in report.problems[:5]:
    print(problem.kind, problem.line, problem.id)
print(report.count("strand_conflict"))
```

问题类型包括 `missing_gene_id`、`missing_transcript_id`、`duplicate_id`、`gene_conflict`、`strand_conflict`、`outside_parent` 等（GFF3 另有 `missing_parent`、`unknown_parent`、`seqid_conflict`、`bad_nesting`、`multiple_parents`）。没有 `gene_id` 的记录、重复的基因 / 转录本记录不进入层级，保存在 `annotation.others` 中。输入不需要排序。

## 表格

注释以列式表格暴露，列是 NumPy 数组或 `StringVector`，行之间用下标互相引用：

```python
genes = annotation.genes
print(len(genes))
print(genes.id[0], genes.name[0], genes.biotype[0])
print(genes.chroms[genes.seqid[0]], genes.start[0], genes.end[0], genes.strand[0])

g = 0
transcripts_of_g = range(genes.transcript_start[g], genes.transcript_end[g])

tx = annotation.transcripts
print(tx.gene[0], tx.length[0], tx.cds_length[0])   # 所属基因行、拼接长度、CDS 总长
print(tx.tags(0), tx.has_tag("MANE_Select")[:10])
features_of_t0 = range(tx.feature_start[0], tx.feature_end[0])
```

| 表 | 内容与排序 | 主要列 |
|---|---|---|
| `genes` | 按 (染色体, start, end) 排序 | `seqid`、`start`、`end`、`strand`、`id`、`name`、`biotype`、`source`、`type`、`transcript_start` / `transcript_end`、`synthesized` |
| `transcripts` | 每个基因的转录本连续排列 | 上述位置列、`gene`、`id`、`name`、`biotype`、`tag_bits`、`tag_values` / `tag_offsets`、`feature_start` / `feature_end`、`length`、`cds_length`、`synthesized` |
| `features` | 每个转录本的子特征连续排列 | `transcript`、`kind`、`start`、`end`、`strand`、`phase`、`type` |
| `exons` | `features` 中的外显子行 | 同上，`feature` 为在 `features` 中的行号 |
| `others` | 不属于层级的记录，按文件顺序 | `type`、`source`、`score`、`id` 等 |

- `strand` 是 `int8`：`1` 为 `+`，`-1` 为 `-`，`0` 为 `.` 或 `?`。
- `seqid` 是 `annotation.chroms` 的下标。
- `features.kind` 是 `FEATURE_KINDS` 的下标：exon、CDS、five_prime_UTR、three_prime_UTR、UTR、start_codon、stop_codon、Selenocysteine、other。
- `phase` 为 0–2，`.` 记为 255。
- `synthesized` 标记由子特征合成的基因或转录本，其范围是子特征的最小 / 最大坐标。
- `tag_bits` 的第 `i` 位对应 `TAGS[i]`（basic、Ensembl_canonical、MANE_Select、CCDS、appris_principal、RefSeq_Select 等），不同来源的拼写变体映射到同一位。

## 查找基因和转录本

```python
gene = annotation.gene("TP53")                # 先按 id 查找，再按名称
print(gene.id, gene.region, gene.strand, gene.transcripts)
print(gene.attributes)                        # 保留的属性

transcript = annotation.transcript("ENST00000269305.9")
print(transcript.gene, transcript.tags, transcript.features)

rows = annotation.gene_rows("TP53", by="name")        # 所有同名基因的行
rows = annotation.transcript_rows("ENST00000269305.9")
attrs = annotation.feature_attributes(int(transcript.features[0]))
```

`gene` / `transcript` 返回第一个匹配；同一 id 在多条染色体上出现（例如 PAR 区基因）时，用 `gene_rows` 取全部。找不到时抛出 `KeyError`。

## 生成查询 key

```python
from omicsweave.query import SequenceReader, TrackReader

mrna_key = annotation.locus("ENST00000269305.9")          # 外显子拼接
cds_key = annotation.locus("ENST00000269305.9", "cds")    # 含终止密码子
utr3_key = annotation.locus("ENST00000269305.9", "utr3")
print(mrna_key)        # Spliced，带转录本链方向，名称为转录本 id

with SequenceReader("genome.fa.gz") as genome, TrackReader("rna.bw") as track:
    mrna = genome.fetch(mrna_key).text          # 5′ → 3′ 的成熟转录本序列
    coverage = track.signal(mrna_key).values    # 沿转录本的信号
```

`locus` 接受转录本 id、名称或行号。`loci(rows, part)` 一次为多行生成 key（默认全部转录本），适合批量构建样本：

```python
selection = annotation.subset(tags=["MANE_Select"])
keys = annotation.loci(selection.transcripts, "cds")
```

## 区域查询

```python
genes_here = annotation.overlaps(("chr17", 7_660_000, 7_700_000), level="gene")
exons_here = annotation.overlaps(("chr17", 7_660_000, 7_700_000), level="exon")
print([annotation.genes.name[int(g)] for g in genes_here])

offsets, rows = annotation.overlaps_batch(
    [("chr1", 0, 1_000_000), ("chr2", 0, 1_000_000)], level="transcript"
)
rows_of_first = rows[offsets[0]:offsets[1]]

nearest, distance = annotation.nearest(("chr17", 8_000_000, 8_000_100))
```

- `level` 为 `"gene"`、`"transcript"` 或 `"exon"`，返回对应表的行号，按起点排序。`exon` 级返回 `annotation.exons` 的行号。
- 未知染色体返回空结果。
- `nearest` 返回所有并列最近的行和距离：重叠为 0，该染色体没有记录为 `-1`。

## 筛选子集

```python
selection = annotation.subset(
    biotypes=["protein_coding"],          # 基因 biotype
    transcript_biotypes=None,             # 转录本 biotype
    tags=["basic", "MANE_Select"],        # 必须同时具有这些已知标签
    tag_values=None,                      # 具有其中任一标签值
    genes=None,                           # 基因 id 或名称
    chroms=["chr1", "chr2"],
)
print(selection.genes, selection.transcripts)
```

各条件同时满足；基因只要有一个转录本满足转录本条件即被选中。`Selection` 可以传给区间推导方法和 `write_bed12`。

## 推导区间

```python
introns = annotation.introns()                           # 每个转录本的内含子，owner 为转录本行
exonic = annotation.exonic_regions()                     # 每个基因外显子的并集，owner 为基因行
promoters = annotation.promoters(upstream=2000, downstream=200)   # 默认按转录本
tss = annotation.tss(level="gene")
tes = annotation.tes()
tes_windows = annotation.tes_windows(upstream=100, downstream=500)
bodies = annotation.gene_bodies(upstream=2000)           # Signac 风格 gene activity 区域
bins = annotation.exon_bins(aggregate=True)              # DEXSeq 式不相交外显子 bin
junctions = annotation.junctions()                       # 已注释的剪接位点
```

| 方法 | 结果 | 默认 |
|---|---|---|
| `introns(subset)` | 相邻外显子之间的区间 | owner：转录本 |
| `exonic_regions(subset)` | 每个基因外显子的并集 | owner：基因 |
| `promoters(upstream, downstream, level, chrom_lengths, subset)` | TSS 上游 `upstream`、下游 `downstream` 的窗口，按链方向 | 2000 / 200，`level="transcript"` |
| `tss(level, subset)` / `tes(level, subset)` | 第一个 / 最后一个转录碱基 | `level="transcript"` |
| `tes_windows(upstream, downstream, ...)` | TES 前 `upstream`（基因体内）、后 `downstream` 的窗口 | 0 / 0 |
| `gene_bodies(upstream, downstream, level, ...)` | 基因体向 TSS 上游、TES 下游扩展 | 2000 / 0，`level="gene"`；ArchR 风格用 `upstream=5000` |
| `exon_bins(aggregate, subset)` | 按覆盖转录本集合切开的不相交外显子 bin | `aggregate=False` |
| `junctions(subset)` | 不重复的内含子，附 donor / acceptor 和所属基因、转录本 | — |

窗口在 0 处截断；提供 `chrom_lengths`（例如 `SequenceReader.chroms`）时也在染色体末端截断。

`IntervalSet` 的列为 `seqid`、`start`、`end`、`strand`、`owner`（无归属为 `2**32 - 1`）和可选 `name`，并带有重叠索引：

```python
rows = promoters.overlaps(("chr1", 10_000, 20_000))               # mode="overlaps"
inside = promoters.overlaps(("chr1", 0, 1_000_000), mode="within")
covering = promoters.overlaps(("chr1", 15_000, 15_001), mode="containing")
owners = [annotation.transcripts.id[int(t)] for t in promoters.owner[rows]]
```

自定义区间（peaks、bins）可用 `IntervalSet.from_arrays(chroms, seqid, start, end, strand, owner, names)` 构建，再交给 `quant.Regions` 计数。

## 染色体名匹配

注释与 BAM / 参考的命名常常不同（`chr1` 与 `1`，`chrM` 与 `MT`）：

```python
from omicsweave.annotation import read_aliases

annotation = Annotation.read("Homo_sapiens.GRCh38.gtf.gz", aliases={"1": "chr1", "MT": "chrM"})
annotation = Annotation.read("genes.gtf", aliases=read_aliases("chromAlias.txt"))

index = annotation.match_chroms(["chr1", "1", "chrUn_xxx"])  # 每个名称在 annotation.chroms 中的下标，-1 为无
```

`match_chroms` 依次尝试同名、给定 aliases 和内置的 `chr1` ↔ `1`、`chrM` ↔ `MT` 规则。`quant` 用同样的规则匹配注释与 BAM。普通序列、信号 reader 不会自动应用这些规则，用注释生成的 key 去查询时，需确保染色体名与数据文件一致（必要时在载入时用 `aliases` 改名）。

## 用于定量

```python
from omicsweave import quant

annotation = Annotation.read("genes.gtf", attributes=None)
task = quant.Task(
    ["reads.bam"], "results",
    annotation=annotation,
    products=[
        quant.Counts.rna_bulk(),
        quant.Counts(quant.Regions(annotation.introns(), feature_type="Introns"), name="introns"),
    ],
)
report = quant.run(task)
```

```sh
ow quant reads.bam -o results --gtf genes.gtf --genes bulk
ow quant cells.bam -o results --gtf genes.gtf --genes single-cell --biotypes protein_coding,lncRNA
```

`GeneExons`（合并外显子，featureCounts / STARsolo `Gene`）、`GeneSpans`（含内含子，STARsolo `GeneFull`）、`Promoters`、`GeneActivity` 都从注释构建特征。`Task.aliases` / `--aliases` 补充内置规则之外的名称映射。详见 [quant](../quant.md)。

## 导出

```python
annotation.write_gtf("copy.gtf", layout="original")
annotation.write_gtf("model.gtf", layout="model")
annotation.write_bed12("transcripts.bed")
annotation.write_bed12("mane.bed", subset=annotation.subset(tags=["MANE_Select"]))
```

| layout | 输出 |
|---|---|
| `original` | 按原文件顺序写出全部记录和注释行；记录由表字段和保留的属性重新序列化，保留属性顺序与引号，属性间距规范化；不保留空行或 CRLF，染色体使用载入后的名称。要求以 `attributes="all"` 从 GTF 载入 |
| `model` | 从表格写出 gene、transcript 和 feature 行 |
| `auto`（默认） | 可以时用 `original`，否则 `model` |

BED12 每行一个转录本。blocks 优先取合并后的 exon 区间；缺少 exon 时使用其他外显子级特征（CDS、UTR、密码子等），仍为空则使用整个转录本范围。thick 区间覆盖 CDS 和起止密码子，并裁剪到 blocks 的整体范围；没有编码特征时 `thickStart = thickEnd = chromEnd`。可再用 `ow convert` 转成 [bigBed](bigbed.md) 按区域查询。这些导出是 `Annotation` 的方法，不属于 `ow convert` 路由。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 用 key 查询时报找不到染色体 | 注释与数据文件命名不同；载入时用 `aliases` 改名 |
| `report` 中有很多 `outside_parent` | 子特征超出了转录本或基因的范围，常见于手工拼接的注释；检查来源文件 |
| 基因数比预期多 | 同一 `gene_id` 出现在多条染色体（如 PAR）时按 (seqid, id) 分别建模 |
| CDS 长度比 `cds_length` 多 3 | `locus(..., "cds")` 包含终止密码子，`cds_length` 是 CDS 行的总长 |
| `write_gtf(layout="original")` 报不支持 | 载入时 `attributes` 不是 `"all"`，或源文件是 GFF3 |
| 内存占用较高 | 不需要原始属性时使用 `attributes=None` |
| `formats.detect` 识别成 BED | GTF 不在 detect 范围内，使用 `Annotation.read` |

## 相关格式

- [GFF3](gff3.md)：用 `ID` / `Parent` 表达层级的注释格式，共享同一套 `Annotation` API。
- [bigBed](bigbed.md)：`write_bed12` 的结果可转换为 bigBed 按区域查询。
- [FASTA](fasta.md)、[bigWig](bigwig.md)：用注释生成的 `Spliced` key 读取转录本序列和信号。

# GFF3

[格式总览](README.md) · [GTF](gtf.md) · [查询](../query.md) · [定量](../quant.md)

GFF3（Generic Feature Format version 3）是描述**基因组特征及其层级**的九列文本格式。它与 [GTF](gtf.md) 的九列结构相同，不同之处在于第九列：GFF3 用 `key=value` 属性，并通过 `ID` 和 `Parent` 显式连接父子特征，因此可以表达任意深度和类型的层级。NCBI RefSeq、Ensembl、众多非模式物种的基因组项目都以 GFF3 发布注释。

在 OmicsWeave 中，GFF3 与 GTF 使用同一个入口 `Annotation.read`，载入后得到同样的“基因 → 转录本 → 特征”列式模型，可以生成 `Spliced` 查询 key、推导内含子 / 启动子等区间，并交给 `quant` 构建基因特征。本文着重说明 GFF3 自身的结构和载入规则；查找、区间推导、定量和导出的完整 API 见 [GTF](gtf.md)。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 基因组特征及其父子关系 |
| 常见后缀 | `.gff3`、`.gff`、`.gff3.gz` |
| 存储 | tab 分隔文本，九列，`##` 指令行 |
| 坐标 | 文件中 1-based、两端闭合；载入后为 0-based 半开 |
| 层级 | `ID` / `Parent` |
| 外部索引 | 不需要；整份载入内存后建索引 |
| 典型来源 | NCBI RefSeq、Ensembl、Phytozome、FlyBase 等，以及 Prokka、BRAKER 等注释流程 |
| OmicsWeave | `Annotation.read` 载入；与 GTF 共享查询、区间推导、quant 特征；可导出为 GTF / BED12 |

## 格式由来

GFF 最早由 Sanger 中心在 1990 年代末提出，第二版（GFF2）衍生出了 GTF。GFF2 的属性语法不统一，层级只能靠约定的属性名表达。Sequence Ontology 项目于 2004 年前后制定 GFF3，统一了属性语法，规定用 `ID` / `Parent` 表达层级，并要求第 3 列使用 Sequence Ontology 术语（`gene`、`mRNA`、`exon`、`CDS`、`ncRNA_gene`……）。

## 文件结构

```text
##gff-version 3
##sequence-region chr1 1 248956422
chr1	RefSeq	region	1	248956422	.	+	.	ID=chr1:1..248956422;Name=1
chr1	BestRefSeq	gene	11874	14409	.	+	.	ID=gene-DDX11L1;Name=DDX11L1;gene_biotype=transcribed_pseudogene
chr1	BestRefSeq	transcript	11874	14409	.	+	.	ID=rna-NR_046018.2;Parent=gene-DDX11L1;Name=NR_046018.2
chr1	BestRefSeq	exon	11874	12227	.	+	.	ID=exon-NR_046018.2-1;Parent=rna-NR_046018.2
chr1	BestRefSeq	exon	12613	12721	.	+	.	ID=exon-NR_046018.2-2;Parent=rna-NR_046018.2
###
```

### 九列

| 列 | 名称 | 含义 |
|---|---|---|
| 1 | seqid | 染色体或 contig 名 |
| 2 | source | 注释来源 |
| 3 | type | Sequence Ontology 类型：`gene`、`mRNA`、`lnc_RNA`、`exon`、`CDS`、`five_prime_UTR` 等 |
| 4 | start | 起点，1-based，包含 |
| 5 | end | 终点，1-based，包含 |
| 6 | score | 分值，通常为 `.` |
| 7 | strand | `+`、`-`、`.`（无方向）、`?`（未知） |
| 8 | phase | CDS 的相位 `0` / `1` / `2`；其他类型为 `.` |
| 9 | attributes | `key=value` 对，以 `;` 分隔 |

### 属性

```text
ID=rna-NM_000546.6;Parent=gene-TP53;Name=NM_000546.6;tag=MANE Select;Note=tumor protein p53%2C transcript variant 1
```

- 以大写字母开头的属性有保留含义：`ID`（文件内唯一标识）、`Parent`（父特征的 ID，可有多个，逗号分隔）、`Name`（显示名）、`Alias`、`Note`、`Dbxref`、`Ontology_term` 等。小写 key 由各数据源自定义，例如 `gene_biotype`、`transcript_id`。
- 多值属性用逗号分隔，例如 `Parent=rna1,rna2` 表示一个外显子被两个转录本共享。
- 保留字符（`;`、`=`、`,`、`&`、tab 等）在值中需要 URL 编码，例如 `,` 写成 `%2C`。
- 一个不连续的特征（例如跨多个外显子的 CDS）可以用多行共享同一个 `ID` 表示。

### 层级

```text
gene              ID=gene-TP53
└── mRNA          ID=rna-NM_000546.6;Parent=gene-TP53
    ├── exon      Parent=rna-NM_000546.6
    ├── CDS       ID=cds-NP_000537.3;Parent=rna-NM_000546.6   （多行共享同一 ID）
    ├── five_prime_UTR
    └── three_prime_UTR
```

转录本这一层的类型名并不固定：`mRNA`、`transcript`、`lnc_RNA`、`ncRNA`、`tRNA`、`V_gene_segment` 等都可能出现。没有转录本的记录（`region`、`repeat_region`、`origin_of_replication` 等）也可以出现在同一文件中。

### 指令行

| 指令 | 含义 |
|---|---|
| `##gff-version 3` | 文件第一行，声明版本 |
| `##sequence-region seqid start end` | 某条序列的范围 |
| `###` | 之前的特征都已完整，后续行不会再引用它们 |
| `##FASTA` | 之后是 FASTA 格式的序列，直到文件结束 |

单个 `#` 开头的行是注释。

### 与 GTF 的对照

| | GTF | GFF3 |
|---|---|---|
| 属性语法 | `key "value";` | `key=value;` |
| 层级 | 共享 `gene_id` / `transcript_id` | `ID` / `Parent` 引用 |
| 层级深度 | 固定为 基因 → 转录本 → 特征 | 任意 |
| 特征类型 | 约定的少数几种 | Sequence Ontology 术语 |
| CDS 与终止密码子 | CDS 通常不含终止密码子，另有 `stop_codon` 行 | CDS 通常包含终止密码子 |
| 附带序列 | 不支持 | 可在 `##FASTA` 之后附带 |

## OmicsWeave 如何载入 GFF3

- **识别**：有 `##gff-version 3` 指令，或第一条记录的属性是 `key=value` 语法时，按 GFF3 解析；也可传入 `dialect="gff3"` 强制指定。纯文本、gzip、BGZF 均可。
- **层级**：沿 `ID` / `Parent` 构建。凡是有外显子级子特征（exon、CDS、UTR 等）的记录都被识别为转录本，无论类型名是什么；其父记录为基因。父记录出现在子记录之后也可以，输入不需要排序。
- **模型 id**：基因优先使用 `gene_id`，转录本优先使用 `transcript_id`；没有对应属性时才使用该记录的 `ID`。`ID` / `Parent` 仍用于连接原始层级，查询和导出则使用模型 id。
- **多个 Parent**：外显子等子特征的多个 Parent 会分别归入各个转录本；转录本有多个父基因时使用第一个，并记录 `multiple_parents` 问题。
- **缺失的记录**：没有自身记录的基因由其转录本合成（`synthesized`），范围是子特征的最小 / 最大坐标。
- **名称与类型**：基因、转录本名称取自 `Name`（`gene_name`、`gene`、`transcript_name` 等也会读取）；biotype 取自 `gene_biotype`、`gene_type`、`biotype`、`transcript_biotype` 等；`tag` 的逗号分隔值作为转录本标签，`MANE Select` 这类 RefSeq 拼写映射到与 GENCODE 相同的已知标签位。
- **其他记录**：不属于基因层级的记录（`region` 等）保存在 `annotation.others` 中。
- **属性值**：按文件原样保存，不做 URL 解码；`%2C` 读出后仍是 `%2C`。
- **指令**：`###` 和其他指令行按注释处理；遇到 `##FASTA` 时停止读取记录，其后的序列被忽略。需要这些序列时，请另外准备 [FASTA](fasta.md) 文件。

| 功能 | 支持情况 |
|---|---|
| 读取 | `Annotation.read`，与 GTF 相同的参数 |
| 查询 key、区间推导、筛选、染色体名匹配 | 与 GTF 完全相同，见 [GTF](gtf.md) |
| quant | `Task.annotation`；CLI 的 `--gtf` 同样接受 GFF3 |
| 导出 | `write_gtf(..., layout="model")` 写成 GTF；`write_bed12` 导出转录本 |
| 写 GFF3 | 不提供 |

GFF3 不通过 `formats.detect` 或区域 reader 打开（`detect` 会把它当作普通文本表），也不需要外部索引。

## 载入与检查

```python
from omicsweave.annotation import Annotation

annotation = Annotation.read("GCF_000001405.40_GRCh38.p14_genomic.gff.gz", attributes=None)
print(annotation.dialect)                      # gff3
print(len(annotation.genes), len(annotation.transcripts))
print(set(annotation.transcripts.type.to_list()))   # mRNA、lnc_RNA、transcript、...

report = annotation.report
print(report.is_clean, report.synthesized_genes)
for problem in report.problems[:5]:
    print(problem.kind, problem.line, problem.id)
```

GFF3 特有的问题类型：

| 问题 | 含义 |
|---|---|
| `missing_parent` | 外显子、CDS 等子特征没有 `Parent` |
| `unknown_parent` | `Parent` 指向的 ID 不存在 |
| `duplicate_id` | 重复 `ID` 的记录中至少一个不是外显子级特征，包括 gene 与 mRNA 共用 `ID`；后者保存在 `others` 中。CDS 等不连续特征可多行共用 `ID` |
| `seqid_conflict` | 子特征与父特征位于不同序列；该记录保存在 `others` 中 |
| `bad_nesting` | 层级不是 基因 → 转录本 → 特征，例如外显子直接挂在已有转录本的基因下、转录本挂在转录本下或 Parent 成环 |
| `multiple_parents` | 转录本有多个父基因，使用第一个 |

此外还有与 GTF 共用的 `strand_conflict`、`outside_parent` 等。`strict=True` 时第一个问题即抛出 `FormatError`。

NCBI 的 GFF3 使用 `NC_000001.11` 这样的 RefSeq 序列名，而参考 FASTA 或 BAM 可能使用 `chr1`。可在载入时改名：

```python
from omicsweave.annotation import read_aliases

annotation = Annotation.read("genomic.gff.gz", aliases=read_aliases("refseq_to_ucsc.txt"))
```

别名文件每行两个名称（原名 新名），`#` 开头的行和空行被跳过，可以由 UCSC `chromAlias` 表取两列得到。

## 查询和定量

```python
key = annotation.locus("NM_000546.6")       # RefSeq 记录含 transcript_id=NM_000546.6 时
cds = annotation.locus("NM_000546.6", "cds")
genes = annotation.overlaps(("NC_000017.11", 7_660_000, 7_700_000), level="gene")
promoters = annotation.promoters(upstream=2000, downstream=200)
```

这里的查询名应与 `annotation.transcripts.id` 中的模型 id 一致：例如 `ID=rna-NM_000546.6;transcript_id=NM_000546.6` 使用 `NM_000546.6`；只有缺少 `transcript_id` 时才使用 `rna-NM_000546.6`。也可按转录本名称查找。`cds` part 拼接该转录本的编码特征（CDS 及起止密码子）；GFF3 的 CDS 通常已包含终止密码子。

```sh
ow quant reads.bam -o results --gtf genes.gff3 --genes bulk
```

CLI 选项名是 `--gtf`，同样接受 GFF3。

## 导出

```python
annotation.write_gtf("genes.gtf")            # GFF3 来源时 auto 即 model 布局
annotation.write_bed12("transcripts.bed")
```

`layout="model"` 从模型写出 GTF 的 gene、transcript 和 feature 行：`gene_id` / `transcript_id` 使用模型 id，即优先读取同名属性、缺省才用原记录的 `ID`；名称和 biotype 写成 GTF 惯用的 key，保留的其他属性附在其后，GFF3 的 `ID` / `Parent` 不写出。`layout="original"` 只适用于从 GTF 载入的注释，不能用于 GFF3 原文回写。当前不提供 GFF3 writer。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 被当作 GTF 解析 | 文件缺少 `##gff-version 3` 且首条记录属性不像 `key=value`；传入 `dialect="gff3"` |
| 大量 `unknown_parent` | 文件被截断或拼接，父记录缺失；检查来源 |
| 属性值中出现 `%2C` 等 | GFF3 的 URL 编码按原样保留，需要时自行 `urllib.parse.unquote` |
| 某些转录本的 type 不是 mRNA | GFF3 允许多种转录本类型；用 `transcripts.type` 或 `subset` 筛选 |
| 染色体名是 `NC_0000xx` | NCBI 命名；用 `aliases` 改为与参考一致的名称 |
| `##FASTA` 之后的序列没有读到 | Annotation 只读取注释记录；序列请使用单独的 FASTA |

## 相关格式

- [GTF](gtf.md)：用 `gene_id` / `transcript_id` 表达层级的注释格式，包含完整的 Annotation API 说明。
- [bigBed](bigbed.md)：`write_bed12` 的结果可转换为 bigBed 按区域查询。
- [FASTA](fasta.md)：读取注释所描述的序列。

# VCF

[格式总览](README.md) · [BCF](bcf.md) · [CSI](csi.md) · [转换](../convert.md) · [查询](../query.md)

VCF（Variant Call Format）是保存**变异位点与样本基因型**的文本格式：SNV、插入缺失、结构变异，以及每个样本在这些位点上的基因型和质量指标。几乎所有变异检测工具（GATK、bcftools、DeepVariant、Strelka……）都输出 VCF，1000 Genomes、gnomAD、ClinVar 等公共资源也以 VCF 发布。

在 OmicsWeave 中，VCF 是**转换源**：它本身不做区域查询，先一次性转换为 [BCF](bcf.md) + [CSI](csi.md)，之后由 `VariantReader` 按区域读取基因型和字段，由 `ConsensusReader` 生成个体化序列，或交给 `quant` 做已知位点屏蔽和等位计数。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 变异位点、等位基因、位点注释（INFO）、每个样本的基因型与指标（FORMAT） |
| 常见后缀 | `.vcf`、`.vcf.gz`（gzip 或 bgzip）、`.vcf.bgz` |
| 存储 | 制表符分隔文本；常用 bgzip 压缩 |
| 坐标 | POS 1-based；OmicsWeave 查询结果统一为 0-based 半开 |
| 随机访问 | 文本本身不支持；bgzip + tabix `.tbi` 可在其他工具中按区域读取 |
| 外部索引 | 在 OmicsWeave 中不使用，转为 BCF 后用 CSI |
| 规范 | hts-specs VCF v4.x（常见 4.2、4.3） |
| 典型来源 | GATK、bcftools call、DeepVariant、公共变异数据库 |
| OmicsWeave | 格式识别、VCF → BCF 转换（可同时建 CSI）；不直接区域查询 |

## 格式由来

VCF 由 1000 Genomes Project 设计（Danecek et al., *Bioinformatics* 2011），目前由 GA4GH 的 hts-specs 维护，与 SAM / BAM、CRAM 属同一套规范体系。它的二进制形式就是 [BCF](bcf.md)。

VCF 的设计重点是“自描述”：header 声明每个字段的含义、类型和长度，数据行只按这些声明书写。这让同一种格式可以承载从单样本 SNV 到数十万样本队列的各种数据，但也意味着读取时必须先解析 header，再逐行解析文本。

## 文件结构

一个 VCF 由三部分组成：`##` 开头的 meta 行、`#CHROM` 开头的列名行、数据行。

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##FILTER=<ID=LowQual,Description="Low quality">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=DB,Number=0,Type=Flag,Description="In dbSNP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	10	rs1	A	G	50	PASS	DP=20;AF=0.5;DB	GT:AD	0/1:10,10	1|1:0,20
chr1	20	.	ACGT	A	20	LowQual	DP=5;AF=0.25	GT:AD	0|1:3,2	./.:.
chr1	30	.	C	T,G	.	.	AF=0.1,0.2	GT:AD	1/2:1,2,3	0:5
```

### Meta 行

| 行 | 作用 |
|---|---|
| `##fileformat=VCFv4.2` | 版本，必须是第一行 |
| `##contig=<ID=...,length=...>` | 声明染色体及长度 |
| `##INFO=<ID,Number,Type,Description>` | 声明位点级字段 |
| `##FORMAT=<ID,Number,Type,Description>` | 声明样本级字段 |
| `##FILTER=<ID,Description>` | 声明 FILTER 列可用的过滤名 |
| `##ALT`、`##reference`、`##source` 等 | 其他说明信息 |

INFO / FORMAT 的 **Number** 说明每条记录有几个值：

| Number | 含义 | 例子 |
|---|---|---|
| 整数 `N` | 固定 N 个值；`0` 表示 Flag | `DP` 为 1 |
| `A` | 每个 ALT 一个值 | `AF` |
| `R` | 每个等位基因（含 REF）一个值 | `AD` |
| `G` | 每种可能基因型一个值 | `PL` |
| `.` | 不定长 | |

**Type** 为 `Integer`、`Float`、`Flag`、`Character`、`String` 之一。

### 数据列

| 列 | 内容 | 说明 |
|---|---|---|
| CHROM | 染色体名 | 应在 `##contig` 中声明 |
| POS | 位置 | 1-based，REF 第一个碱基的位置 |
| ID | 变异 ID | 多个用 `;` 分隔，`.` 表示无 |
| REF | 参考等位 | 参考基因组上的碱基 |
| ALT | 替代等位 | 多个用 `,` 分隔 |
| QUAL | 质量分 | Phred 标度，`.` 表示缺失 |
| FILTER | 过滤状态 | `PASS`、过滤名列表或 `.` |
| INFO | 位点字段 | `KEY=VALUE` 用 `;` 分隔，Flag 只写 KEY |
| FORMAT | 样本字段的键 | 如 `GT:AD:GQ`，规定后续样本列的顺序 |
| 样本列 | 每个样本的值 | 按 FORMAT 顺序，用 `:` 分隔 |

### 基因型 GT

GT 用等位编号表示：`0` 是 REF，`1` 是第一个 ALT，`2` 是第二个 ALT，依此类推。

| 写法 | 含义 |
|---|---|
| `0/1` | 杂合，未定相（不知道两个等位分别在哪条染色体上） |
| `0\|1` | 杂合，已定相：第一条单倍型为 REF，第二条为 ALT |
| `1/1` | 纯合 ALT |
| `1/2` | 两个不同的 ALT，多等位位点上常见 |
| `./.` | 缺失 |
| `0` | 单倍体（如男性 X、线粒体） |

### 插入缺失与锚定碱基

VCF 中 REF 和 ALT 必须非空，因此插入缺失会带上事件**前面一个参考碱基**作为锚定碱基：

```text
chr1	20	.	ACGT	A	...   # 删除 CGT：REF 覆盖 1-based 20–23
chr1	50	.	T	TGA	...   # 在 50 之后插入 GA
```

在 OmicsWeave 的查询结果中，第一条记录为 `start=19, end=23`（0-based 半开，覆盖整个 REF）。因此，一个只与 `CGT` 重叠的查询区域也会返回这条删除。

### 结构变异：symbolic ALT 与 END

大片段变异不写出完整序列，而用尖括号的符号等位，并用 INFO/END 或 SVLEN 给出跨度：

```text
chr1	100	sv1	N	<DEL>	60	PASS	SVTYPE=DEL;END=200	GT	0/1
```

常见符号有 `<DEL>`、`<DUP>`、`<INV>`、`<INS>`、`<CNV>`，gVCF 中还有表示“非变异区块”的 `<*>` / `<NON_REF>`。断点（breakend）用 `N[chr2:123[` 这类写法表示。

### 压缩与 tabix

`.vcf.gz` 可能是普通 gzip，也可能是 bgzip（BGZF）。只有 bgzip 压缩并配套 tabix `.tbi` 时，htslib 类工具才能按区域读取文本 VCF。扩展名不能区分两者。

## 关键概念

### 为什么查询前要转 BCF

文本 VCF 的每一行都要切分字符串、按 header 解释类型；多样本队列中一行可能有数十万个样本列。[BCF](bcf.md) 把同样的信息预先编码成类型化的二进制，查询时可以直接拷贝数值，只解码需要的字段（例如只要 GT）。OmicsWeave 只对 BCF 做区域查询，VCF 统一先转换一次。

### Header 必须完整

BCF 把字段名编码为 header 字典中的编号，因此 VCF 中用到的 contig、FILTER、INFO、FORMAT 都必须在 header 中声明。很多 VCF 缺少 `##contig` 行（如 ClinVar），可以用 `add_contigs` 补声明；其他缺失的声明需先修正源文件。

### 排序

建 CSI 索引要求每个 contig 的记录连续，且同一 contig 内 POS 不减；contig 块的出现顺序不必与 header 中的声明顺序一致，但离开某个 contig 后不能再次返回。OmicsWeave 的 VCF → BCF 转换保留输入顺序，不排序；当前也没有 VCF 排序路由，未排序的 VCF 应先用外部工具（如 `bcftools sort`）排序。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 格式识别 | 以 `##fileformat=VCF` 开头的文本识别为 `Vcf`；gzip / bgzip 压缩识别为 `GzipVcf` |
| 输入压缩 | 纯文本、普通 gzip、BGZF 均可作为转换输入 |
| 区域查询 | 不支持；`query.open("x.vcf")` 会报错，并提示对应的 `ow convert` 命令 |
| VCF → BCF | 支持，输出与 `bcftools view -Ob` 相同的 BCF 记录流，按 htslib 的方式切 BGZF 块 |
| 同时建索引 | `index=True` / `--index` 写 `<bcf>.csi`（`min_shift` 固定为 14，与 `bcftools index` 默认相同） |
| 补 contig | `add_contigs=True` / `--add-contigs`：额外扫描一遍 CHROM 列，按首次出现顺序补 `##contig=<ID=...>` |
| 压缩等级 | `level` / `--level`，0–12 |
| 版本 | 转换要求首行以 `##fileformat=VCF` 开头，缺失时报错；通过此前缀检查后，版本号无法解析时才按 VCFv4.2 处理 |
| tabix `.tbi` | 不读取，也不需要 |

限制：

- 转换不排序；`index=True` 遇到未排序记录会报错。
- 未声明的 contig 默认报错（`contig 'chr2' is not declared in the header`）；未声明的 INFO / FORMAT / FILTER 需先修正源文件，`add_contigs` 只处理 contig。
- `add_contigs` 补出的 contig 没有长度，`VariantReader.chroms` 中对应值为 `None`。
- 转换不修改记录：不做左对齐、多等位拆分、规范化或去重。
- 当前不提供 BCF → VCF 导出，也没有从内存写 VCF / BCF 的 writer。

## 使用

### 转换为 BCF 并查询

```sh
ow convert cohort.vcf.gz cohort.bcf --index
```

```python
from omicsweave import convert
from omicsweave.query import VariantReader

report = convert.run("cohort.vcf.gz", "cohort.bcf", index=True, threads=4)
print(report.target, report.companions)  # Bcf('cohort.bcf') (Csi('cohort.bcf.csi'),)
print(report.records)                    # 转换的变异记录数

with VariantReader("cohort.bcf") as variants:
    print(variants.samples)
    gt = variants.genotypes(("chr1", 0, 100_000))
    print(gt.gt.shape)  # (记录数, 样本数, 最大倍性)
```

查询结果的字段、形状和缺失值见 [BCF](bcf.md)。同一个 BCF 还可用于 consensus、quant 的已知变异屏蔽和等位计数。

### header 缺少 contig

```sh
ow convert clinvar.vcf.gz clinvar.bcf --add-contigs --index
```

`--add-contigs` 会多读一遍输入，只补 contig 声明，不排序。

### 未排序的 VCF

```sh
# 不建索引时可以转换，但得到的 BCF 不能按区域查询
ow convert unsorted.vcf.gz unsorted.bcf

# 先用外部工具排序，再转换并建索引
bcftools sort unsorted.vcf.gz -Oz -o sorted.vcf.gz
ow convert sorted.vcf.gz sorted.bcf --index
```

### 显式指定格式

内容识别通常足够。文件名不规范、或想避免识别开销时，可以显式给出格式对象：

```python
from omicsweave.formats import GzipVcf

convert.run(GzipVcf("calls.txt.gz"), "calls.bcf", index=True)
```

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| `VariantReader` 拒绝 `.vcf.gz` | VCF 文本不做区域查询；先 `ow convert x.vcf.gz x.bcf --index` |
| `contig ... is not declared in the header` | header 缺 `##contig`；加 `--add-contigs` |
| `INFO 'XX' is not declared in the header` | 转换遇到 header 未声明的 INFO 字段；先修正源 VCF |
| 建索引报 `unsorted positions` | 输入未按坐标排序；先排序，或不建索引 |
| 已有 `.tbi` 却仍需转换 | OmicsWeave 的变异查询基于 BCF + CSI，不使用 tabix |
| 查询区域返回了起点在区域之前的记录 | 记录与区域按 REF 跨度（或 END）判断重叠，删除和 SV 会跨入区域 |

## 相关格式

- [BCF](bcf.md)：VCF 的二进制形式，OmicsWeave 中变异查询的格式。
- [CSI](csi.md)：BCF 的区域索引。
- [FASTA](fasta.md)、[2bit](2bit.md)：consensus 所需的参考序列。

# 定量分析

[按区域查询](query.md) · [格式转换](convert.md)

`quant` 从按坐标排序的 BAM 生成覆盖度、计数矩阵和实验相关统计。一次任务可以同时请求多个产物，共用对每个输入 BAM 的一次扫描。

例如，同时生成 RNA-seq 的覆盖度和基因计数：

```sh
ow quant sample.bam -o results \
  --coverage --genes bulk --gtf genes.gtf \
  --matrix mtx -t 4 -m 2G
```

这会生成 `results/coverage.bw`、`results/genes/` 下的计数矩阵与 `counts.tsv`，以及任务统计 `results/summary.json`。文中命令都需要替换为自己的输入路径。

## 工作方式

```text
坐标排序 BAM                      ┌─► Coverage     ─► coverage.bw
(+ 注释 / 参考 / BCF)             ├─► Counts       ─► genes/matrix.*、barcodes、features
        │                         ├─► Junctions    ─► junctions/junctions.tab
        ▼                         ├─► Methylation  ─► methylation/*.tsv.gz、*.bw
  一次顺序扫描，每条记录 ─────────┼─► Contacts     ─► hic/contacts.pairs.gz、contacts.hic
  分发给所有请求的产物            └─► ...
                                         ▼
                                  summary.json（所有产物的统计）
```

- **一次扫描，多个产物。** BAM 只解压、解码一次，每条记录交给所有请求的产物；产物各自过滤、各自计数，互不影响。同时请求覆盖度、基因计数和剪接位点，和单独请求其中一个的扫描成本几乎相同。
- **从已比对数据开始。** 输入是比对器产出的坐标排序 BAM，条形码、UMI 等信息来自 BAM 标签。quant 不做比对、条形码纠错或细胞判定。
- **内存与输出解耦。** 每个产物在缓冲预算内累加，超出时写入临时文件，结果与是否 spill、线程数无关。
- **输出标准格式。** 轨道是 bigWig，矩阵是 CellPack / Matrix Market / CSV，接触是 4DN pairs 和 .hic，片段与胞嘧啶表是 bgzip + tabix。下游工具可以直接读取，也可以用 OmicsWeave 的 reader 按区域查询。

## 按实验类型速查

| 实验 | 典型命令 | 需要的额外输入 |
|---|---|---|
| 任意测序的覆盖度轨道 | `ow quant in.bam -o out --coverage` | — |
| 链特异 RNA-seq 覆盖度 | `ow quant in.bam -o out --coverage --strands split --library 2 --scale cpm` | — |
| bulk RNA-seq 基因计数 | `ow quant a.bam b.bam -o out --gtf genes.gtf --genes bulk --pairs --library 2` | GTF / GFF3 |
| 10x 风格 scRNA-seq | `ow quant cells.bam -o out --gtf genes.gtf --genes single-cell --matrix mtx` | GTF / GFF3；BAM 含 `CB` / `UB` |
| RNA velocity | `... --genes single-cell --velocity` | GTF / GFF3 |
| 剪接位点 | `ow quant in.bam -o out --junctions --gtf genes.gtf --reference genome.fa.gz` | 参考（motif）、注释（是否已注释） |
| scATAC-seq | `ow quant atac.bam -o out --atac --peaks peaks.bb --gtf genes.gtf` | 可选 peaks bigBed、注释 |
| ChIP-seq / CUT&Tag | `ow quant chip.bam -o out --chip-peaks peaks.bb` 或 `--chip-bins 5000` | peaks bigBed |
| SLAM-seq 等代谢标记 | `... --labeling --labeled 2 --conversion 'T>C' --reference genome.fa.gz` | 参考或 `MD` 标签；可选已知变异 BCF |
| WGBS / EM-seq / ONT 甲基化 | `ow quant methyl.bam -o out --methylation --reference genome.fa.gz` | 参考（参考比对方式需要） |
| Hi-C / Micro-C | `ow quant hic.bam -o out --hic --hic-resolutions 5000,25000` | BAM 含 `MC` / `MQ` |
| 等位计数 / ASE | `ow quant in.bam -o out --alleles sites.bcf --ase sample1 --gtf genes.gtf` | 位点 BCF + CSI；ASE 需定相 |

每一项的计数定义、默认过滤和输出文件见下文对应小节。

## 支持的分析

| 目的 | Python 产物 | CLI 入口 | 主要输出 |
|---|---|---|---|
| 覆盖度 | `Coverage` | `--coverage` | 一个或分链的 bigWig |
| 单细胞 / bulk RNA 计数 | `Counts.rna_single_cell` / `rna_bulk` | `--genes single-cell` / `bulk` | 基因计数矩阵，bulk 另有 `counts.tsv` |
| 任意区域计数 | `Counts` | 常用场景见 ATAC、ChIP 选项 | region、peak、bin 等特征的矩阵 |
| scATAC | `Atac` | `--atac` | fragments、索引、细胞指标、peak / tile / gene activity 矩阵 |
| ChIP-seq / CUT&Tag | `Counts.chip` | `--chip-peaks` / `--chip-bins` | 片段计数矩阵和表 |
| 剪接位点 | `Junctions` | `--junctions` / `--junction-matrix` | STAR 九列 junction 表，可选矩阵 |
| RNA velocity 计数层 | `Counts(..., velocity=True)` | `--genes ... --velocity` | `spliced`、`unspliced`、`ambiguous` |
| 代谢标记 | `Labeling`、`Counts(..., labeled=Labeled(...))` | `--labeling`、`--genes ... --labeled N` | 转换直方图、错配统计、labeled 计数层 |
| DNA 甲基化 | `Methylation` | `--methylation` | 胞嘧啶表、beta / coverage 轨道、M-bias、区域矩阵 |
| Hi-C / Micro-C | `Contacts` | `--hic` | 4DN pairs、原始计数 `.hic`、可选细胞指标 |
| SNV 等位计数与 ASE | `Alleles` | `--alleles`、`--ase` | `ref` / `alt` / `other`，基因 `hap1` / `hap2` / `both` |

这些功能从已有比对开始。细胞条形码和 UMI 使用 BAM 中保存的值；不做条形码纠错、UMI 纠错、细胞判定、从头 peak calling 或从头变异检测。velocity 产物是计数层，代谢标记产物是观测统计，均不替代下游模型拟合。

## 准备输入

- 输入为一个或多个按坐标排序的 [BAM](formats/bam.md)。quant 顺序扫描 BAM，不要求 `.bai` / `.csi`。
- SAM、CRAM 先[转换](convert.md)成 BAM；quant 不直接接收它们。
- 剪接表和 Hi-C 接触合并多个 BAM 时，参考序列名称、长度和排列顺序必须一致。需要每个样本一行时使用 `rows="@file"`；以条形码为行时，同名条形码会合并，合并样本前应先区分条形码。
- `Coverage` 和 `Atac` 每次只支持一个 BAM；多个样本应分别运行。计数矩阵等产物可以接收多个 BAM。
- 基因、启动子、gene activity、velocity、基因 ASE 等分析需要 [GTF](formats/gtf.md) / [GFF3](formats/gff3.md) 注释。
- 参考序列接受普通 FASTA + `.fai`、bgzip FASTA + `.fai` + `.gzi`，或 2bit。代谢标记和参考比对甲基化可在没有参考文件时使用读段的 `MD` 标签；剪接 motif 需要参考文件。
- 已知变异屏蔽和等位计数使用 BCF + `.csi`。

未排序 BAM 先处理为：

```sh
ow convert unsorted.bam sorted.bam --sort -t 4 -m 2G
```

CLI 的 `--gtf` 同时接受 GTF 与 GFF3。注释和 BAM 的 `chr1` / `1`、`chrM` / `MT` 可以匹配；额外名称映射用 `--aliases aliases.tsv`，每行是注释名称和 BAM 名称。

## Python：组合一次任务

```python
from omicsweave import quant
from omicsweave.annotation import Annotation

annotation = Annotation.read("genes.gtf", attributes=None)
task = quant.Task(
    inputs=["sample.bam"],
    out_dir="results",
    annotation=annotation,
    matrix="mtx",
    products=[
        quant.Coverage(),
        quant.Counts.rna_bulk(unit="pair"),
        quant.Junctions(matrix=quant.JunctionMatrix.bulk()),
    ],
)
report = quant.run(task, threads=4, memory="2G")
print(report.files)       # 每个文件的 kind 和 path
print(report.records)     # 扫描到的有参考位置的记录数，不是分子数
print(report.summary)     # 与 summary.json 相同的统计
```

`Task.products` 中的产物各自有过滤和计数规则；加入覆盖度不会改变基因计数的规则。需要多个同类型产物时，通过 `name=` 为它们指定不同输出名称。

`Task` 的其他字段：

| 字段 | 作用 |
|---|---|
| `annotation` | 基因类产物、ATAC gene activity 与 TSS 指标使用的 `Annotation` |
| `aliases` | 注释名 → BAM 名的额外染色体映射；`chr1` ↔ `1`、`chrM` ↔ `MT` 已内置 |
| `matrix` | 矩阵存储：`"cpk"`（默认）、`"mtx"`、`"csv"` |
| `reference` | 比对所用参考（FASTA + 索引或 2bit），用于代谢标记、剪接 motif 与甲基化调用 |
| `variants` | 已知位点 BCF（带 CSI），计数转换和甲基化时屏蔽这些位点 |

## 输出目录与 summary.json

以开头的 RNA-seq 命令为例，输出目录是：

```text
results/
├── coverage.bw              # Coverage
├── genes/                   # Counts.rna_bulk
│   ├── matrix.mtx.gz        # 主计数层（--matrix mtx）
│   ├── barcodes.tsv.gz      # 观测行：这里是每个输入 BAM
│   ├── features.tsv.gz      # 特征：id, name, type, chrom, start, end
│   └── counts.tsv           # featureCounts 风格的表
└── summary.json             # 本次任务的全部统计
```

每个产物写在以产物名命名的位置（覆盖度是 `<name>.bw`，矩阵类是 `<name>/` 目录），用 `name=` 可以改名，从而在一个任务中放多个同类产物。

`summary.json` 记录了每个输入和每个产物的统计，适合在流程中自动检查：

```json
{
  "omicsweave": "0.1.0",
  "inputs": [{"path": "sample.bam", "placed_records": 2}],
  "products": [
    {
      "product": "counts",
      "unit": "read",
      "rows": 1, "features": 1, "nonzero": 1,
      "assigned": 2,
      "unassigned_no_features": 0,
      "unassigned_ambiguity": 0,
      "filtered_flag": 0, "filtered_mapq": 0, "filtered_multimapping": 0,
      "files": ["genes/barcodes.tsv.gz", "genes/features.tsv.gz", "genes/matrix.mtx.gz", "genes/counts.tsv"]
    }
  ]
}
```

（上例省略了部分字段。）`filtered_*` 统计被过滤条件排除的记录，`unassigned_*` 按原因统计没能分配到特征的单位，与 featureCounts 的 summary 类似。Python 中 `report.summary` 是同一份内容，`report.files` 列出每个输出文件的 `kind`（如 `bigwig`、`mtx`、`barcodes`、`features`、`fragments`、`tbi`、`json`）和路径。

## 计数矩阵的含义和存储

矩阵的两个维度是观测单位和特征：观测单位可以是条形码、read group 或输入文件；特征可以是基因、peak、bin 或自定义区间。

| `matrix=` / `--matrix` | 文件 | 存储方向 |
|---|---|---|
| `cpk`，默认 | `matrix.cpk` | 条形码 / 样本 × 特征，CellPack `uint32` |
| `mtx` | `matrix.mtx.gz` | 特征 × 条形码 / 样本，Matrix Market |
| `csv` | `matrix.csv.gz` | 特征为行，条形码 / 样本为列，稠密表 |

目录同时包含 `barcodes.tsv.gz` 和 `features.tsv.gz`。后者的六列为 `id, name, type, chrom, start, end`，坐标是 0-based 半开。Matrix Market 的矩阵下标从 1 开始；`counts.tsv` 为 featureCounts 风格，`Start` 从 1 开始。

附加层使用层名，例如 `spliced.mtx.gz`、`labeled.mtx.gz`；甲基化使用 `methylated.*` / `total.*`，等位计数使用 `ref.*` / `alt.*` / `other.*`。这些产物不一定有名为 `matrix.*` 的主层。

需要对接常见 10x 分析工具时，可以选择 `matrix="mtx"`。`csv` 会生成稠密矩阵，大型单细胞数据应优先选择稀疏存储。

## 覆盖度

```sh
ow quant reads.bam -o coverage --coverage
ow quant reads.bam -o coverage-stranded --coverage --strands split --library 2 --scale cpm
```

```python
product = quant.Coverage(
    name="rna",
    unit="bases",
    strands="split",
    strand="reverse",
    scale="cpm",
)
```

| `unit` | 一个读段或片段贡献什么 |
|---|---|
| `bases`，默认 | CIGAR `M` / `=` / `X` 覆盖的碱基；`deletions=True` 时也计 `D`，不计 `N` |
| `span` | 比对的整个参考跨度 |
| `five-prime` / `three-prime` | 读段的 5′ / 3′ 端位置 |
| `fragment` | 一对读段对应的片段跨度，从最左 mate 计一次；单端读段按 span |

默认不缩放。`scale="cpm"` 按计入的读段或片段数缩放到每百万；也可给数值因子。合链输出 `<name>.bw`，分链输出 `<name>.forward.bw` 和 `<name>.reverse.bw`。

默认过滤排除 unmapped、secondary、QC-failed、duplicate（`0x704`），不限制 MAPQ，不要求 `NH=1`。重叠 mates 的碱基可以重复贡献，未提供 `samtools depth -s` 式重叠去除；fragment 过滤以最左 mate 为准。

## RNA 与区域计数

```sh
ow quant cells.bam -o scrna --gtf genes.gtf --genes single-cell --matrix mtx
ow quant a.bam b.bam -o bulk --gtf genes.gtf --genes bulk --pairs --library 2
```

| 预设 | 观测行 | 特征与计数方式 | 默认规则 |
|---|---|---|---|
| `Counts.rna_single_cell()` | `CB` | 合并的基因外显子，按 `UB` 精确去重 | forward，`NH=1`，排除 `0x704` |
| `Counts.rna_bulk()` | 每个 BAM | 基因外显子，按 read；可用 `unit="pair"` | unstranded，`NH=1`，仅排除 unmapped，另写 `counts.tsv` |
| `Counts.atac(features)` | `CB` | Tn5 切点，按片段位置去重 | 命中多个特征时均计入 |
| `Counts.chip(features)` | 每个 BAM | 片段 | 排除 `0x704`，命中多个特征时均计入，另写表 |

预设可以通过关键字覆盖。是否属于链特异文库要按实际实验设置，不应只凭测序类型推断。

### 选择特征

| 特征对象 | 含义 | 需要的额外输入 |
|---|---|---|
| `GeneExons()` | 每个基因合并后的外显子 | 注释 |
| `GeneSpans()` | 整个基因跨度，包含内含子 | 注释 |
| `Promoters(upstream=2000, downstream=200)` | 按基因链方向的 TSS 窗口 | 注释 |
| `GeneActivity(upstream=2000, downstream=0)` | 基因体加上下游扩展 | 注释 |
| `Peaks("peaks.bb")` | bigBed 中每条区间 | bigBed，BED 需先转换 |
| `Bins(5000)` | 每条参考序列上的固定大小窗口 | BAM 的参考序列信息 |
| `Regions(intervals)` | 一个 `IntervalSet` 的每条区间 | 自定义区间或注释生成的区间集 |

例如，对注释中的内含子计数：

```python
introns = annotation.introns()
product = quant.Counts(
    quant.Regions(introns, feature_type="Introns"),
    name="introns",
    rows="@file",
    unit="read",
)
```

### 选择单位、去重与分配

`unit` 支持 `read`、`pair`、`fragment`、`tn5-fragment`、`cut-site`。read pair 按两端的比对块分配；fragment 按片段跨度分配；Tn5 片段做端点校正；cut-site 对校正后的两个切点分别分配。

`dedup="none"` 不去重，`"umi"` 在每个观测行、每个特征内按 UMI 精确值去重，`"position"` 按片段端点去重。UMI 不作编辑距离聚类；同一 UMI 在不同基因上可以分别计数。大写 ACGT 且长度不超过 31 的 UMI 使用无碰撞编码，其他 UMI 使用 63-bit 指纹，未额外消解碰撞；按 read name 标识的分子使用 64-bit 指纹。

一般以至少 1 bp 重叠分配。默认重叠多个特征的单位视为 ambiguous；`multi_overlap=True` 时计入每个命中特征。pair 模式优先选择两端共同支持的特征；两端分别过滤，一端被过滤时另一端可按单端计数。

`rows` 可设为任意标签名（如 `"CB"`）、`"@RG"` 或 `"@file"`。`whitelist` 是每行一个条形码的文件，既限制条形码，也确定输出顺序，包含没有读段的指定条形码。

未提供 fractional 计数、最大重叠分配或读段延伸；与其他工具比较时应统一这些定义以及过滤条件。

## scATAC 与 ChIP / CUT&Tag

```sh
ow quant atac.bam -o atac-results --atac --peaks peaks.bb --gtf genes.gtf --matrix mtx
ow quant chip.bam -o chip-results --chip-peaks peaks.bb --coverage
```

`Atac()` 默认要求 proper pair、MAPQ ≥ 30，排除 unmapped、mate-unmapped、secondary、QC-failed、supplementary。duplicate 标记的记录保留，再按条形码与片段端点折叠。过滤只检查用于生成片段的最左 mate，不检查另一端的 MAPQ。

片段使用 `[pos + 4, pos + |TLEN| - 5)`；两个切点为 `start` 和 `end - 1`。同一条形码、同一染色体、相同端点的片段合并，`fragments.tsv.gz` 的第五列是支持它的 read pair 数。

默认输出：

- `atac/fragments.tsv.gz` 和 `.tbi`；`--no-fragments` 可关闭。
- `atac/barcode_metrics.csv`：片段、peak 内切点、TSS 及侧翼切点等指标。
- `atac/tiles/`：默认 5000 bp tile 的切点矩阵；`--tiles 0` 可关闭。
- 指定 `--peaks` 时的 `atac/peaks/` 切点矩阵。
- 提供注释时的 `atac/gene_activity/` 片段矩阵，默认基因体加上游 2000 bp；`--gene-activity none` 可关闭。

ChIP / CUT&Tag 的 `--chip-peaks` 和 `--chip-bins N` 分别对 peaks 和固定大小 bin 计片段，不做 peak calling。

## 剪接位点与 velocity

```sh
ow quant rna.bam -o splice --junctions --gtf genes.gtf --reference genome.fa.gz
ow quant cells.bam -o velocity --gtf genes.gtf --genes single-cell \
  --velocity --junction-matrix single-cell --matrix mtx
```

`Junctions` 从 CIGAR 的 `N` 提取内含子，输出 `junctions/junctions.tab`，九列依次为染色体、内含子首末碱基（1-based，包含两端）、链、motif、是否已注释、unique 计数、multi-mapping 计数、最大 overhang。默认一对读段跨同一位点计一次；`--single-mates` 改为每端分别计数。

有参考时识别 GT/AG 等 motif，没有参考时 motif 为 0。有注释时判断位点是否已注释。输出全部观测位点，不做 STAR `--outSJfilter*` 的计数、overhang 或内含子长度过滤。`JunctionMatrix.single_cell()` 生成按 `CB`、UMI 去重的矩阵；`bulk()` 按输入文件计数。

velocity 为基因计数增加 `spliced`、`unspliced`、`ambiguous` 三层。外显子和内含子证据按分子合并后分类；主计数层仍使用配置的特征分配规则。实现未提供 velocyto 的两遍 validated-intron 流程，不保证与其他 velocity 工具整套流程完全等价。

## 代谢标记

```sh
ow quant slam.bam -o slam --gtf genes.gtf --genes single-cell \
  --labeling --labeled 2 --reference genome.fa.gz \
  --variants known.bcf --conversion 'T>C' --matrix mtx
```

`Labeling` 输出 `labeling/conversions.tsv.gz` 和 `labeling/mismatches.tsv`，分别记录每条读段的源碱基覆盖数 / 转换数直方图，以及按 mate、测序 cycle、参考碱基分类的错配统计。总体替换比例写入 `summary.json`。

`Labeled(min_conversions=2)` 为基因矩阵添加 `labeled` 层；主层 `matrix.*` 仍是 total。转换按转录方向定义：T>C 对正链是基因组 T>C，对负链是 A>G。默认碱基质量阈值为 20，可屏蔽 `Task.variants` 的已知位点。

参考来自 `Task.reference`，否则使用读段的 `MD`。read pair 合并两端转换位置，同一位置计一次；UMI 分子只要有一条读段达到阈值即标记，不把所有读段的转换数累加判定。单独的 `Labeling` 直方图始终按每条读段统计。

不做读段两端 trimming，也不拟合新旧 RNA 比例或转换率模型。`--conversion 'T>C'` 必须加引号，避免 shell 将 `>` 解释为重定向。

## DNA 甲基化

```sh
ow quant methyl.bam -o methyl-results --methylation \
  --reference genome.fa.gz --contexts CpG,CHG,CHH --merge-cpg
```

| `source` / `--methylation-source` | 调用来源 |
|---|---|
| `auto`，默认 | 逐记录依次尝试 Bismark 标签、MM/ML、参考比对 |
| `bismark` | `XM` / `XG`，按标签给出的上下文和状态 |
| `modifications` | MM/ML 中的 5mC `C+m` 概率 |
| `reference` | 与参考或 MD 重建的参考比较，判断 C>T / G>A |

默认只统计 CpG，最低碱基质量为 5，MM/ML 概率阈值为 0.5；重叠 mates 丢弃 read 2 在 read 1 区间内的调用。`merge_cpg=True` 合并 CpG 两条链，CHG / CHH 不作链合并。

输出包括 `methylation/cytosines.tsv.gz` 与 `.tbi`、`beta.bw`、`coverage.bw`、`mbias.tsv`。胞嘧啶表七列为 `chrom, start, end, context, strand, methylated, unmethylated`。beta 是甲基化比例，coverage 是总调用数；没有覆盖的位置不列出。M-bias 仅报告，不自动截断或过滤读段。

区域矩阵用两层保存计数，适合下游自行计算比例：

```python
product = quant.Methylation(
    table=False,
    bigwig=False,
    matrices=[quant.MethylMatrix("bins", quant.Bins(100_000), rows="@file")],
)
```

每胞嘧啶表和 bigWig 只支持单个 BAM；多个 BAM 时关闭 `table` 和 `bigwig`，或在 CLI 使用 `--no-cytosines`。区域矩阵输出 `methylated.*` 与 `total.*`，支持 bins、peaks、基因、启动子和自定义区域。MM/ML 当前仅输出 5mC，不单独输出 5hmC 等其他修饰，也不自动估计概率阈值。

## Hi-C / Micro-C

```sh
ow quant hic.bam -o hic-results --hic --hic-resolutions 5000,25000 --genome hg38
```

输入需要 mate 标签 `MC` 和 `MQ`。如果标签缺失，应先按 samtools fixmate 所需的名称分组流程补齐，再恢复为坐标排序 BAM。

接触取自两个 mate 的 5′ 位置，默认两端 MAPQ ≥ 30，排除非主比对、QC-failed 和 duplicate 等记录，再按端点、链及条形码精确去重。`min_distance` 过滤距离过近的 cis 接触，`self_ligation` 对 inward / outward 接触设置距离阈值。

输出 `hic/contacts.pairs.gz` 和 `hic/contacts.hic`。`--hic-resolutions none` 可只输出 pairs，`--no-pairs` 可只输出 hic。`--hic-cells --barcode-tag CB` 另写 `hic/cell_metrics.csv`。

只使用主比对，不处理 walks 或嵌合多段解析；pairs 的 pair type 为 `UU`，不生成 pairix 索引。`.hic` 只有原始计数，不生成 expected、归一化向量或 FRAG 分辨率。它与 pairtools 完整解析及模糊坐标去重流程不同。

## 等位计数与 ASE

```sh
ow quant cells.bam -o allele-results --alleles sites.bcf --matrix mtx
ow quant bulk.bam -o ase-results --alleles sites.bcf --allele-mode bulk \
  --ase sample1 --gtf genes.gtf --library 0
```

位点输入为 BCF + `.csi`，仅统计 REF / ALT 都是单个 ACGT 碱基的二等位 SNV。多等位、indel 等记录跳过并记入统计，不从 BAM 发现新位点。

默认单细胞按 `CB` 和 `UB` 计分子，bulk 按文件和 read name 计片段；默认 MAPQ ≥ 20、碱基质量 ≥ 20。输出 `alleles/ref.*`、`alt.*`、`other.*` 以及 `variants.tsv.gz`。所有支持的位点保留为特征，不做最小计数或最小等位频率筛选。

一个分子在同一位点出现冲突观测时，会在出现的每个等位层各计一次，不做多数表决。因此 `ref + alt + other` 不一定等于独立分子总数。

`AlleleGenes(sample="sample1")` / `--ase sample1` 使用该样本的定相杂合位点，将基因单位分为 `hap1`、`hap2`、`both`，输出到 `alleles/genes/`。需要注释和已有定相，不进行单倍型推断。

## 过滤、链方向与资源

Python 用每个产物的 `filter=quant.Filter(...)` 设置 `require`、`exclude`、`min_mapq` 和 `unique`。`unique=True` 要求 `NH=1`；没有 NH 时按 1 处理。不同预设的默认过滤不同，替换整个 Filter 时应明确给出需要保留的条件。

链方向 `unstranded`、`forward`、`reverse` 分别对应 CLI `--library 0 / 1 / 2`：forward 表示 read 1 与转录本同向，reverse 表示相反；read 2 按其配对关系判断。CLI 共享选项会应用到每个具有对应选项的产物；需要对产物分别设置时使用 Python。

先查看计划而不扫描全部读段：

```sh
ow quant sample.bam -o results --coverage --dry-run -t 4 -m 2G
```

```python
plan = quant.plan(task, threads=4, memory="2G", tmp_dir="scratch")
print(plan.inputs, plan.budget, plan.memory)
report = plan.run()
```

计划读取 BAM header 并准备特征和所需辅助数据，不能代替对全部记录的检查。默认使用进程可用的 CPU 核心，缓冲预算为 4 GiB。超过各产物的缓冲份额时使用磁盘临时文件；`tmp_dir` 默认是输出目录，应预留可用空间。预算不是 RSS 硬上限。

任务完成后通过 `report.files` 和 `summary.json` 检查输出文件、过滤、未分配原因与跳过记录。输出先写临时文件，整项计算成功后发布；多个文件的发布不是文件系统事务，应在 `run()` 成功返回后使用整组结果。

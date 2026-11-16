# BAM 定量分析

[文档首页](README.md) · [格式转换](conversion.md) · [基因注释](annotation.md)

## 一次扫描生成多个产物

```sh
# 单细胞 RNA：UMI 基因矩阵和覆盖度
ow quant possorted.bam -o out --gtf genes.gtf --genes single-cell --coverage -t 4 -m 4G
# bulk RNA：按 read pair 计数，反链文库
ow quant a.bam b.bam -o counts --gtf genes.gtf --genes bulk --pairs --library 2
# ATAC：fragments、peak / tile / gene-activity 矩阵
ow quant atac.bam -o atac --gtf genes.gtf --atac --peaks peaks.bb --tiles 5000
# 剪接和 velocity
ow quant possorted.bam -o splice --gtf genes.gtf --reference genome.fa.gz \
    --genes single-cell --velocity --junctions --junction-matrix single-cell
```

输入必须是带 `SO:coordinate` header 且实际按坐标排序的 BAM。量化按序扫描，不要求 BAI / CSI；CRAM、SAM 需先[转换](conversion.md)。产物的过滤和计数单位由各自选项决定，不应把所有矩阵解释成同一种“read count”。

```python
from omicsweave import quant
from omicsweave.annotation import Annotation

task = quant.Task(
    inputs=["possorted.bam"],
    out_dir="out",
    products=[quant.Coverage(), quant.Counts.rna_single_cell()],
    annotation=Annotation.read("genes.gtf", attributes=None, threads=4),
    matrix="mtx",
)
plan = quant.plan(task, threads=4, memory="4G")
print(plan.memory, plan.needs, plan.inputs)
report = plan.run()  # 或 quant.run(task, ...)
print(report.files, report.summary["products"])
```

`plan` 读取 header 并建立特征空间；注释、白名单、变异面板也可能占用常驻内存，因此它不是完全无 I/O 的校验。`Plan.run()` 直接运行这份计划，不再重新规划。默认内存预算 `"4G"`，默认使用进程可用的全部核心。产物都是 frozen dataclass，默认值和预设（`Counts.rna_single_cell` / `rna_bulk` / `atac` / `chip`、`JunctionMatrix.single_cell` / `bulk`、`Alleles.single_cell` / `bulk`）来自 Rust 层，规划时由 Rust 逐项校验。矩阵默认 `cpk`；`mtx` 和 `csv` 可显式选择。`--dry-run` / `quant.plan` 不生成定量产物。

| 产物 | Python 配置 | 主要输出 |
|---|---|---|
| 覆盖度 | `Coverage` | bigWig，可分链和缩放 |
| 特征计数 | `Counts`（特征：`GeneExons`、`GeneSpans`、`Promoters`、`GeneActivity`、`Peaks`、`Bins`、`Regions`；`subset=Annotation.subset(...)` 限定基因与转录本） | 基因、区域或 bin 矩阵，bulk 另有 `counts.tsv` |
| ATAC | `Atac` | BGZF fragments + `.tbi`、指标、区域矩阵 |
| 剪接 | `Junctions` / `JunctionMatrix` | 位点表及可选矩阵 |
| velocity / 标记层 | `Counts.rna_single_cell(velocity=True, labeled=Labeled(...))` | spliced / unspliced / ambiguous / labeled |
| 代谢标记 | `Labeling` | 转换直方图和 mismatch 表 |
| 甲基化 | `Methylation` / `MethylMatrix` | 胞嘧啶表、beta / coverage、M-bias、矩阵 |
| Hi-C | `Contacts` | pairs、hic、细胞指标 |
| 等位基因 | `Alleles` / `AlleleGenes` | SNV ref / alt / other 与基因 hap1 / hap2 / both |

普通 FASTA 需要 `.fai`，bgzip FASTA 需要 `.fai` / `.gzi`，2bit 无需索引。已知变异和 SNV 位点输入为 BCF + CSI。没有参考时，部分调用可从 `MD` 标签恢复参考；无法恢复的读段会计入相应的未分配统计。Hi-C 使用 `MC` / `MQ` mate 标签；准备后须再次保证坐标排序。全部 CLI 选项见 `ow quant --help`。

## 输出和兼容性

计数产物各有 `barcodes.tsv.gz`、`features.tsv.gz` 及矩阵。CPK 按 barcode 存行；Matrix Market 为 **feature × barcode**，两者方向不同。其他层与主矩阵共用 feature / barcode 表。`summary.json` 包含每个产物的计数和未分配原因，解释结果时应一起检查。

UMI 按精确值去重，不纠错、不做细胞判定。Velocity 不实现 velocyto 的两遍 validated-intron 模型；甲基化不做自动 M-bias trimming；Hi-C 不处理 supplementary walk，也不输出归一化向量；SNV 模块不做变异发现或单倍型定相。下文保留每个产物的完整计数定义和与参照工具的差异。

## 执行模型与计数定义

`omicsweave-quant` 把已比对、按坐标排序的 BAM 一次扫描成覆盖度轨道、计数矩阵（含 velocity 与标记层）、单细胞 ATAC fragments、剪接位点表、代谢标记统计、DNA 甲基化、Hi-C 接触和 SNV 等位基因计数。边界：不做比对，不链接、不复用 STAR；与格式 crate 和 `ow convert` 不同，这一层做分析（过滤、特征分配、计数、去重），格式转换的职责边界不限制定量分析。

### 一次扫描，多个状态机

每个产物是一个累加器（`Accumulator`），所有累加器共享同一次扫描：

1. `source`：打开时检查 header 的 `@HD SO:coordinate`，扫描时检查 (tid, pos) 单调不减；任一不满足即报错，并提示 `ow convert in.bam sorted.bam --sort`。BGZF 块在线程池上解压（`omicsweave_bam::Scanner`，建在 `omicsweave_bgzf::Scanner` 上），读取线程把连续记录切成批（与 convert 相同，64 KiB 到 4 MiB，随预算变化），遇到第一条未放置的记录（tid = -1，按序排在最后）即结束。多个 BAM 依次扫描，各自是一个独立的行来源。BAI / CSI 只会用于限定区域，目前不读。
2. `prelude`：`ReadView` 在第一次被问到时才计算并缓存派生量，同一条记录的多个累加器共享：FLAG / MAPQ / `NH` 过滤（`Filter`）、可选字段查找、文库链特异性、比对块（`M`/`=`/`X`，可选并入 `D`）、内含子（`N`）、5′ / 3′ 端、片段（由最左 mate 的 TLEN 得出）、Tn5 校正片段和切点、比对段（`segments()`：每个 `M`/`=`/`X` 的参考起点、读段起点、长度，供 `query_at` / `reference_at` 双向查位置）、读段下的参考碱基（`reference_bases()`：比对区间两侧各多取 2 bp，供序列上下文）、读段区间内的已知变异位置（`masked()`），以及对齐碱基（`aligned()`：每个 `M`/`=`/`X` 位置的参考位置、读段位置、读段碱基、参考碱基、质量、是否落在已知变异上）。参考碱基来自任务的参考基因组（`reference`：bgzip / 普通 FASTA 加 `.fai` / `.gzi`，或 2bit，经 `omicsweave_fasta`），没有参考或该染色体不在参考里时用 `MD` 标签；已知变异来自 BCF + CSI（`omicsweave_vcf::BcfReader`）。参考和变异在计划时按名字（含别名与 `chr1` ↔ `1`）对到每个输入的染色体，长度不符即报错；每批一个只向前滑动的窗口（按记录位置取序列、查询变异，前后预取 256 kb），每个碱基每批约取一次。mate 配对分两级：批内按 (read name, 本端位置, mate 位置) 配对（`Pairs`），配不上的半对在合并时交给 `Mates`；当合并前沿越过 mate 应在的位置仍未出现时，该半对按单端计。
3. `assign`：特征空间（`FeatureSpace`）是一个 `omicsweave_gtf::IntervalSet` 加上「行 → 列」映射：基因的合并外显子、基因全长、启动子 / gene activity 窗口（`Derive` 按 BAM 的染色体长度截断）、bigBed 的全部 peak（BED 不直接读，先用 `ow convert` 转成 bigBed）、固定大小的 bin。任务内按构造方式去重注册；每批每个空间一个 `SweepCursor`，染色体名用 `match_seqids` / `match_names` 对应（`chr1` ↔ `1`、别名表）。一条记录对同一空间、同一查询（块、片段、切点等）、同一链要求只分配一次。
4. `accum`：见下文。
5. `model` 与 `sink`：产物是与文件格式无关的数据（`Track`、`CountMatrix`、`FragmentTable`、`BarcodeMetrics`、`JunctionTable`、`Histogram`、`MismatchTable`、`Cytosine` 流、`Mbias`、`ContactTable`、`VariantTable`；矩阵的列是 `FeatureList`，特征表或按需生成行的位点列表），大的产物是只读一次的有序流（记录可能在 spill 文件里），由 sink 写成文件。

### Accumulator

每个累加器有两种状态：批状态（`Batch`）在 worker 上对一批连续记录 `observe`，结束时 `seal`（排序、折叠、批内配对，仍在 worker 上）；总状态（`Total<'p>`）在单个合并线程上按批次顺序 `merge`，合并时得到该批最后一条记录的坐标键（frontier），据此判断哪些结果已经最终确定。`finish` 把总状态变成产物，`deliver` 交给 sink 并返回 `summary.json` 中的一段。类型擦除的 `Stage` / `Merging` 让任务持有任意组合。总状态在 `RunCtx` 给出的预算份额内自己 spill（外排序 `Sorter`、有序缓冲 `Spool`）。`Needs`（cigar、tags、mates、assign、sequence、quality、reference、variants）声明累加器读什么，计划时取并集；只有并集含 reference / variants 时批才建参考 / 变异窗口，计划时为每个在途批多计 2 MiB / 1 MiB。

新增一个产物只需写一个累加器：声明 `Needs`、`estimate` 内存（不可 spill 的固定部分、每批在途部分、最小缓冲和缓冲权重）、批如何观察一条记录、批如何合并、交付什么。

### 预算与 spill

复用 convert 的 `Budget { threads, memory, tmp_dir }` 和 `MemoryPlan`：`plan` 只读 BAM header 并构建特征空间，给出 `fixed`（特征表、索引、白名单、注释）、`pipeline`（在途的块和批）、`buffer`（各累加器的排序缓冲和 spool）。`fixed + pipeline + 各累加器最小缓冲` 超出预算时直接报错；其余缓冲按权重分给各累加器。总状态超出份额时 spill 到 `tmp_dir`（默认输出目录），不能 spill 的部分（无白名单时的条形码字典、等待配对的 mate）超出份额时报错并提示给白名单或增加内存。所有输出先写到临时文件，全部产物成功后才改名；产物在各自线程上并行写出，一起发布。

### 确定性

批次按序合并，所有折叠满足交换律与结合律（整数求和、取最大、去重计数）；行号按条形码首次出现的记录顺序编号（与批的切分无关）；覆盖度的 run 跨批合并成极大段。所以输出与线程数、批大小、是否 spill 都无关，测试逐字节检查（`summary.json` 只记录与此无关的统计量）。

### 覆盖度（`Coverage`）

单位：比对碱基（`samtools depth` 的语义，可选计入 `D`，即 `-J`；`N` 不计）、整段 span、5′ 端、3′ 端、片段（每对只从最左 mate 计一次，`[pos, pos + |TLEN|)`，可设最大长度；单端读段按 span，非最左 mate 不计）。链：合并，或按文库类型拆成正负两条（read 2 视作 read 1 的链，反向文库再翻转）。过滤默认排除 `0x704`，可设最小 MAPQ 和只计 `NH == 1`。缩放：原值、CPM（每百万计入的读段或片段）或给定因子。批把区间端点记成排序的事件；总状态保留前沿之后的事件，把之前的扫成等值 run 写入 spool，内存只随前沿处的深度增长，与基因组大小无关；最后写成 bigWig（`BigWigWriter`，带 zoom，可选 sidecar）。

### 计数矩阵（`Counts`）

一个累加器，维度可组合：

- 特征（列）：上面的任一特征空间。
- 单位：读段；read pair（两端比对块的并，featureCounts `-p --countReadPairs`；不带 `multi_overlap` 时按 featureCounts 的投票：被两端都覆盖的特征优先，并列才算模糊）；片段（最左 mate 的 TLEN 区间，可做 Tn5 校正）；Tn5 切点（校正片段的首、末碱基，各自分配）；分子 = 读段 + UMI 去重。
- 行：标签（`CB` 或任意标签）、read group（各 header 的 `@RG` 依次合并）、输入文件；可给白名单（行按白名单顺序，含没有读段的条形码）。
- 去重：无；UMI（每个 (行, 特征) 中不同的 `UB` 精确值计一次，跨批全局去重：条目带 UMI 进入外排序，读出时折叠）；位置（同一行、同一片段位置只计一次，ATAC 的重复片段）。
- 分配规则同 featureCounts：至少 1 bp 重叠；链特异性 0 / 1 / 2；命中多个列的单位记为模糊（`multi_overlap` 时计入每一列）；多重比对（`NH > 1`）默认不计。

批输出单元格条目 `(row, col, layer, sub, value)`（`sub` 是 UMI 或片段位置），在批内排序折叠；总状态把条目送进外排序，矩阵按行读出时再按 `Fold` 折叠：计数层和标记层按 `Reduce`（求和、按 sub 去重计数、按 sub 取最大），velocity 证据层按位或后分类到三层；内存由排序缓冲界定，与非零元数量无关。

预设：单细胞 RNA（genes × UMI × `CB`，STARsolo `Gene` 的原始矩阵：全部条形码、不做细胞判定、`CB` / `UB` 按原样使用、正链文库、只计唯一比对、排除 `0x704`）；bulk RNA（genes × read 或 pair × 文件，featureCounts 默认：除未比对外全部记录都计、多重比对按 `NH` 排除，另写 featureCounts 格式的 `counts.tsv`）；scATAC（peaks / tiles × 切点 × `CB`，gene activity（基因体加上游 2 kb，Signac）× 片段 × `CB`，均按位置去重）；ChIP / CUT&Tag bulk（bins 或 peaks × 片段 × 文件，片段计入它重叠的每个特征）。

### 单细胞 ATAC（`Atac`）

定义（10x Genomics `fragments.tsv`、ArchR / Signac）：

- 片段：通过过滤（默认要求 `0x1 | 0x2`，排除未比对、mate 未比对、secondary、QC fail、supplementary，MAPQ ≥ 30；重复保留）、两端比对到同一染色体的 read pair，取最左 mate：`[pos + 4, pos + |TLEN| - 5)`，即最左 mate 5′ 端 +4 到另一 mate 5′ 端 -5。
- 切点：片段的首碱基和末碱基（`start`、`end - 1`）。
- 重复：同一条形码、同一 (染色体, start, end) 是一个片段，第五列是其 read pair 数。
- 每条形码指标（`barcode_metrics.csv`）：read pair 数、片段数、有切点落在 peak 内的片段数、peak 内切点数、TSS ±50 bp（101 bp）内切点数、距 TSS 1900-2000 bp 的两个 100 bp 侧翼内切点数；ArchR 的 TSS enrichment 为 `(tss_cuts / 101) / (tss_flank_cuts / 200)`。

片段在批内排序折叠，总状态把前沿之前的片段（位置不可能再有重复的）按序写入 spool，最后写成 BGZF 的 `fragments.tsv.gz`，同时用 `BinningIndex`（BAI 参数）建索引，换成 tabix 的头部（`tbx_conf_bed`、按出现顺序的序列名）得到 `.tbi`。ATAC 的各矩阵与片段文件使用相同的过滤和片段定义，因此各目录的 `barcodes.tsv.gz` 一致。

### 剪接（`Junctions`）

定义（STAR `SJ.out.tab`）：

- 剪接位点是 `N` 操作的参考区间。过滤默认排除未比对、QC fail、supplementary（STAR 的嵌合片段），重复和 secondary 保留；多重比对的每条比对都计。
- 一对 read 同时跨过同一位点只计一次（`pairs`，按 read pair 计数同样的配对；配不上的半对单独计）。
- unique / multi：`NH` 为 1 / 大于 1 的跨越读段（对）数。
- overhang：读段在该位点两侧相邻比对段中较短的一段，比对段是不被插入、缺失、clip、内含子打断的 `M`/`=`/`X`；表中取最大值。
- motif：参考序列上内含子首两个与末两个碱基：1 `GT/AG`、2 `CT/AC`、3 `GC/AG`、4 `CT/GC`、5 `AT/AC`、6 `GT/AT`，其余（及无参考时）为 0；链：奇数 motif 为 `+`，偶数为 `-`，非经典位点有注释时取注释的链，否则 0。annotated：等于注释中某个转录本的内含子（`Annotation::junctions()`，`find_any`）。

批输出 `(tid, start, end, unique, multi, overhang)` 条目，批内排序折叠，总状态进外排序；结束时折叠成 `JunctionTable`（按位置排序，motif 在合并线程上逐位点取 4 个碱基），写成 `junctions/junctions.tab`（STAR 的九列）。可选的 junction × 行矩阵（`JunctionMatrix`：行键、白名单、UMI 去重、额外过滤；STARsolo `SJ` 预设为 `CB` × UMI、只计唯一比对）：列是表中全部位点（按表顺序，含只有多重比对读段的位点），条目先按位点排序（列号在所有位点已知后才确定），读出时编号再进第二个外排序按行读出；与表写在同一个 `junctions/` 目录下。多个 BAM 要求参考序列相同且顺序一致。

### Velocity（`Counts::velocity`）

基因矩阵上的 spliced / unspliced / ambiguous 三层（velocyto / STARsolo `Velocyto` 规则）：读段的比对块与恰好一个基因的全长（含内含子，按文库链）重叠才归属该基因，重叠多个基因不计；对该基因每个包含读段全部比对块的转录本，每个块都落在同一外显子内为 exonic，否则（块在内含子里或跨外显子-内含子边界）为 intronic；没有转录本包含读段则不计。读段的证据位（exonic、intronic）按分子（有 UMI 时 `(行, 基因, UMI)`；无 UMI 时单个读段，read pair 按 read name 哈希）在 `value` 中按位或（`Fold`：该层按位或，其余层用 `Reduce`），读出时只有 exonic 为 spliced、只有 intronic 为 unspliced、两者都有为 ambiguous。与计数层的基因分配相互独立（计数层仍按外显子）。

### 代谢标记（`Labeling`、`Counts::labeled`）

SLAM-seq / scSLAM-seq / TimeLapse / 4sU：

- 读段的转录链按文库：read 2 视作 read 1 的链，反向文库翻转（`Labeling` 把非链特异文库按正链读；基因矩阵在非链特异文库下用所计基因的链）。转换按转录本方向给出（默认 `T>C`，可配置任意 `X>Y`）：正链转录本在基因组上是 `T>C`，负链是 `A>G`。
- 只计对齐碱基中参考与读段碱基都是 `A`/`C`/`G`/`T`、质量不低于 `min_quality`（默认 20；没有质量的记录全部通过）、不在已知变异上的位置；覆盖度是转录本方向为源碱基（`T`）的位置数，转换数是其中读成目标碱基（`C`）的位置数。`convert` 一次算出两条链的结果（转换位置列表），供各累加器共用。
- 基因矩阵的 `labeled` 层：计数层计入的单位（读段、read pair、UMI 分子）中转换数不少于 `min_conversions` 的；read pair 合并两端的转换位置（同一位置计一次），分子只要有一条读段达到阈值就算标记。条目与计数层同 `sub`、同 `Reduce`，在单独的层号上，所以 `matrix.*` 就是 total。
- `Labeling` 产物：每读段 `(T 覆盖度, T>C 数)` 直方图 `labeling/conversions.tsv.gz`（稀疏三列）；`labeling/mismatches.tsv`：按 mate、测序 cycle（从读段测序 5′ 端起，含 hard clip）、转录本方向参考碱基，统计四种读段碱基；`summary.json` 给出总转换率和十二种替换的总体比率。没有参考也没有 `MD` 的读段计入 `unassigned_no_reference`。

### DNA 甲基化（`Methylation`）

WGBS、EM-seq、单细胞亚硫酸氢盐、带碱基修饰标签的长读段（Bismark / MethylDackel / modkit 的定义）：

- 调用（call）是读段上的一个胞嘧啶：参考位置（`C`，反链胞嘧啶取 `G` 的位置）、链、上下文、状态。来源按记录选择（`CallSource::Auto`）：有 `XM` 用 Bismark 的 `XM`（`Z/z` CpG、`X/x` CHG、`H/h` CHH，大写为甲基化，`U/u` 与 `.` 不是调用；链取 `XG`：`CT` 正链，`GA` 反链）；否则有 `MM`/`ML` 用碱基修饰（`C+m` 列表里每个 C 的概率 `(ML + 0.5) / 256` 不小于 `min_probability` 为甲基化、不大于 `1 - min_probability` 为未甲基化，中间不调用；未列出的 C 在默认 / `.` 模式下为未甲基化，`?` 模式下不调用；`MN` 与读段长度不符或有 hard clip 而无 `MN` 时视为标签不合；链取读段的）；否则对参考比对（bwa-meth / MethylDackel）：原始正链读段（`YD:Z:f`、`XG:Z:CT`，都没有时 read 1 正向 / read 2 反向）在参考 `C` 上读到 `C` 为甲基化、`T` 为未甲基化，原始负链读段在参考 `G` 上读到 `G` / `A`，用 `aligned()`。
- 上下文取参考（有参考基因组用基因组，否则用 `MD` 推出的、同一读段相邻对齐位置的参考碱基；`MM` 调用没有参考时用读段自身碱基）：其后为 `G` 是 CpG，否则 `H` 后接 `G` 是 CHG，否则 `HH` 是 CHH（反链在互补链上向左看）；后两位未知的胞嘧啶不调用。只调用质量不低于 `min_quality`（默认 5；没有质量的记录全部通过）、不在任务已知变异上的碱基。
- 重叠 mate 只计一次（Bismark `--no_overlap`）：同一染色体、相距小于 64 kb 的一对在 `Pairs` / `Mates` 中配对（两端判断相同，配不上的半对过期后按单端计），read 2 落在 read 1 比对区间内的调用丢弃。
- M-bias：按上下文、mate、测序方向的读段位置（1-based，含 hard clip）统计甲基化 / 未甲基化调用，在丢弃重叠之前；写 `mbias.tsv`，并写入 `summary.json`。
- 每胞嘧啶表与轨道：批内按 (tid, pos, code) 折叠，总状态保留前沿之后的胞嘧啶，并且不越过等待配对的半对中最小的调用位置（这些半对的第一个调用键放在一个计数多重集中），之前的写入 spool；最后一次读出同时写 `cytosines.tsv.gz`（chrom、start、end、context、strand、methylated、unmethylated，0-based 半开，BGZF，tabix 按 BED 建索引，与 ATAC fragments 共用 `IndexedBed`）、`beta.bw` 和 `coverage.bw`；`merge_cpg` 把一个 CpG 两条链合成一行 `[C, C + 2)`、链为 `.`。只支持单个 BAM（关掉表和轨道后矩阵可以跨多个 BAM，如每细胞一个 BAM）。
- 区域矩阵（`MethylMatrix`，allcools / scBS-map 风格）：特征为任一特征空间（基因、启动子、gene activity、peak、bin），行键同计数矩阵；每个调用按位置分配（`assign_at`，同一批的扫描游标），条目 `(row, col, layer, 0, 1)` 以 `Fold::direct(Sum)` 折叠成 `methylated` / `total` 两层（各层文件都按层命名：`methylated.cpk`、`total.cpk`）。

### Hi-C / Micro-C（`Contacts`）

输入是按坐标排序、带 mate 标签（`MC`、`MQ`，`samtools fixmate -m`）的 BAM（4DN pairs、pairtools 的定义）：

- 接触取自一对中坐标序在前的 mate（同位置取 read 1），用它的 RNEXT/PNEXT、mate 链（0x20）、`MC`（反向 mate 的 5′ 端需要它的 CIGAR）和 `MQ`，不缓存 mate。两端是两个 mate 的 5′ 位置（1-based），按染色体序号再按位置排成 `(chr1, pos1) <= (chr2, pos2)`。
- 过滤：默认要求 paired，排除未比对、mate 未比对、secondary、QC fail、supplementary、duplicate；两端 MAPQ 都不低于 `min_mapq`（默认 30，mate 的取 `MQ`，没有 `MQ` 只在阈值为 0 时通过）；同一染色体上距离小于 `min_distance` 的丢弃；`self_ligation` 时距离小于它的 inward（`+ -`，dangling end）和 outward（`- +`，self-circle）接触丢弃，同向接触不可能是自连接产物，保留。
- 去重：两端位置、链（和条形码）相同的接触只留文件序第一个。
- 批输出变长的 `Contact`（含 read name）进外排序，键 (chr1, chr2, pos1, pos2, strand1, strand2, 行, 输入, 记录号)，所以重复相邻且顺序与批切分无关；读出时去重，同时写 `contacts.pairs.gz`（BGZF，头部含 `#chromsize`，按 BAM 参考序）、把每个接触按每个分辨率分 bin（与 Juicer `pre` / `HicBuilder` 相同：位置截到染色体长度、`bin = pos / res`、染色体内取上三角）送进各分辨率自己的外排序（键 (c1, c2, block, y, x)，`ABSORBS` 折叠计数），最后粗到细逐张 map、按 block 序经 `HicWriter` 的 `begin_matrix` / `add_block` / `end_matrix` 写出 `contacts.hic`（v9，原始计数，可选 sidecar）。convert 的 pixel → `.hic` 阶段不是公开 API，所以这里直接调用 hic writer，业务逻辑（过滤、去重、两端定义）都在 quant。
- 单细胞：按 `barcode_tag` 统计每个条形码的 pairs、duplicates、contacts、cis、trans、cis ≥ 1 / 10 / 20 kb（pairtools stats），写 `cell_metrics.csv`，同时写入 `summary.json` 的 `cells`。

### 等位基因计数（`Alleles`）

cellsnp-lite mode 1 / vartrix：

- 位点是 BCF（带 CSI）中的二等位 SNV（REF、ALT 都是单个 `ACGT`），其他记录跳过；计划时整体读入（每个位点 12 字节，`SitePanel`），按文件序全部作为列（列表行按需生成，不为每个位点存字符串：`SiteFeatures` 实现 `FeatureList`）。
- 读段在位点上有质量不低于 `min_quality`（默认 20）的对齐碱基即观测到它：等于 REF 为 `ref`，等于 ALT 为 `alt`，其他 `ACGT` 为 `other`，`N` 不计；位点落在缺失或内含子中不计。
- 单位：UMI（`umi_tag` 精确值，没有 UMI 的读段不计）或 read name（一对的两端只计一次），即 `(行, 位点, 等位, sub)` 上的 `Reduce::Distinct`；一个单位的读段不一致时在它出现的每个等位上各计一次。行键同计数矩阵。
- 每位点表 `variants.tsv.gz`：矩阵写出时按列累加各层（`write_with` 的 `column_sums`），给出 REF / ALT / other 单位数和 `alt / (ref + alt)`。
- 基因水平 ASE（`AlleleGenes`）：取一个样本（默认第一个）定相杂合（`0|1`、`1|0`）位点，读段读到单倍型 1 的等位为第一位证据、另一等位为第二位；读段的比对块与恰好一个基因的外显子重叠（按文库链）时，证据按单位在 `(行, 基因, sub)` 上按位或，读出时分成 `hap1`、`hap2`、`both` 三层（与 velocity 共用两位证据的折叠）。

### 输出布局

输出目录下：`<name>.bw`（或 `.forward.bw` / `.reverse.bw`）；计数产物 `<name>/`：`barcodes.tsv.gz`、`features.tsv.gz`（id、name、type、chrom、start、end 六列，Cell Ranger ARC 的列）、`matrix.cpk`（CellPack，`u32`，每行一条记录，按行顺序追加；默认）或 `matrix.mtx.gz`（Matrix Market，特征 × 条形码，同 Cell Ranger）或 `matrix.csv.gz`（稠密，一行一个特征），其他层为 `<layer>.cpk` 等（`spliced`、`unspliced`、`ambiguous`、`labeled`），bulk 另有 `counts.tsv`；`atac/fragments.tsv.gz`（+ `.tbi`）、`atac/barcode_metrics.csv`、`atac/peaks/`、`atac/tiles/`、`atac/gene_activity/`；`junctions/junctions.tab`（有矩阵时同目录为 10x 目录，特征类型 `Splice Junction`，id 为 `chrom:start-end`，0-based 半开）；`labeling/conversions.tsv.gz`、`labeling/mismatches.tsv`；`methylation/cytosines.tsv.gz`（+ `.tbi`）、`methylation/beta.bw`、`methylation/coverage.bw`、`methylation/mbias.tsv`、`methylation/<矩阵>/`（`methylated.*`、`total.*`）；`hic/contacts.pairs.gz`、`hic/contacts.hic`、`hic/cell_metrics.csv`；`alleles/`（`ref.*`、`alt.*`、`other.*`）、`alleles/variants.tsv.gz`、`alleles/genes/`（`hap1.*`、`hap2.*`、`both.*`）；`summary.json`（输入、`Needs`、每个产物的计数与各类未分配原因、spill 字节数、文件列表）。压缩文本都是 BGZF（与 gzip 兼容），在线程池上压缩。CellPack 通过 git 依赖（固定 rev）接入，不修改其代码。

### 对齐测试


- 覆盖度：`refs/samtools` 中 vendor 的 `bam2depth.c`（`samtools depth -aa`，含 `-J`、`-Q`），逐碱基比较；其余单位、链拆分对照测试中按定义直接计算的逐碱基计数；1 线程与 8 线程、spill 与不 spill 的 bigWig 逐字节一致。
- bulk 计数：`refs/subread`（Subread 2.2.1 的 featureCounts，GPL-3.0，`publish = false`，只作测试依赖，不进入 wheel）在随机合成的 BAM 和 GTF（重叠基因、两条链、剪接读段、多重比对、重复、secondary / supplementary）上逐基因比较，单端 `-s 0/1/2`、`-O`，以及 `-p --countReadPairs`；`Assigned`、`Unassigned_NoFeatures`、`Unassigned_Ambiguity`、`Unassigned_MultiMapping`、`Unassigned_Unmapped` 也一致。
- 单细胞、read group、白名单、bins、三种矩阵格式：对照按定义直接计算的结果。
- ATAC：片段文件逐行、peak / tile 切点矩阵、gene activity、指标对照按上述定义的直接计算；`.tbi` 与 htslib `tbx_index_build` 对同一文件建出的索引结构一致（bin 顺序除外：htslib 按哈希表顺序写）。
- 剪接：随机合成的读段（成对、多重比对、重复、带插入缺失 clip 的剪接读段，另加正好剪在注释内含子上、两端都跨同一位点的 read pair），基因组随机生成后在大部分位点上写入六种经典 motif；`junctions.tab` 逐行对照按上述定义的直接计算（参考由 htslib `bgzf` + `fai_build3` 写成 bgzip FASTA，另写 2bit，两者输出逐字节一致；无参考时对照无 motif 的版本），junction × `CB` 矩阵和条形码顺序、velocity 三层（单细胞 UMI 与 bulk 读段）同样对照直接计算；1 线程与 7 线程、spill 与不 spill 的文件逐字节一致。
- 代谢标记：从随机基因组按 CIGAR（插入、缺失、内含子、soft / hard clip）抽取读段，正负链转录本、单端与 read pair，注入 `T>C`（转录本方向）、测序错误、`N`、样本 SNP 和随机质量；直方图、mismatch 表、单细胞（UMI，阈值 2）与 bulk（非链特异，按基因链）的 `labeled` 层对照直接计算，带与不带 SNP 掩码（BCF 由 htslib 写出并用 `bcf_index_build3` 建 CSI）；bgzip FASTA、2bit、`MD` 标签（测试按基因组生成）三种参考来源及不同线程数的输出逐字节一致。

- 甲基化：从生成的基因组抽取定向文库的 read pair（两条原始链、每个分子逐个胞嘧啶随机甲基化、测序错误、`N`、插入缺失与 clip、重叠 mate、被过滤的记录、已知 SNP），参考调用、`MD` 调用、`XM` 调用的 `cytosines.tsv.gz`（含 CpG 合并）、beta / coverage bigWig、M-bias、bin × `CB` 与 bin × 文件矩阵对照按上述定义的直接计算，`.tbi` 与 htslib `tbx_index_build` 一致；`MM`/`ML` 长读段（5hmC 列表在前、显式 / 隐式模式、两个概率阈值、有无参考）同样对照直接计算；1 线程与 6 线程、spill 与不 spill 的文件逐字节一致。
- Hi-C：生成的 read pair（cis / trans、两条链、clip 与插入缺失、`MC` / `MQ` 偶尔缺失、任一端低 MAPQ、带标志的记录、未标记的重复、条形码）；`contacts.pairs.gz` 与每条形码统计对照直接计算，`contacts.hic` 与 `HicBuilder` 对同一组接触写出的文件逐字节一致，并用 `HicReader` 读回 25 kb 的每张 map 对照直接分 bin；距离过滤、不去重、MAPQ 0 的 bulk 设置同样；不同线程数、spill 的文件逐字节一致。
- 等位基因：生成的 read pair（定相样本两条单倍型的分子、PCR 拷贝、错误、`N`、低质量、位点上的内含子与缺失、被过滤的记录），BCF 由 htslib 写出（含要跳过的插入和多等位记录）；单细胞 UMI 与 bulk read pair 的 ref / alt / other 层、`variants.tsv.gz`、hap1 / hap2 / both 基因层对照直接计算；不同线程数、spill 的文件逐字节一致。

### 与参照工具的有意偏差


- featureCounts：基因按 (seqid, `gene_id`) 区分（同一 id 在两条染色体上是两个基因）；`counts.tsv` 每个基因一组坐标而不是逐外显子的分号列表；pair 模式下两端各自过滤（featureCounts 对一对整体判断 MAPQ、`NH`），一端被过滤时另一端按单端计；不提供 fractional 计数、`--largestOverlap`、`--fracOverlap`、读段延伸。
- STARsolo / Cell Ranger：UMI 不做纠错，也不过滤映射到多个基因的 UMI：读段各自分配，一个 UMI 的读段落在不同基因时在每个基因各计一次；不做细胞判定。
- 10x fragments：没有 `#` 头部注释行；过滤只看最左 mate（另一 mate 的 MAPQ 不检查）。
- 覆盖度 fragment 单位同样只看最左 mate；不提供 `samtools depth -s` 的重叠去除。
- STAR `SJ.out.tab`：不做 `--outSJfilter*` 过滤（overhang、计数、与其他位点的距离、内含子长度），输出全部位点；非经典位点的链取注释的链；没有参考时 motif 全为 0（不读 `jM` / `XS` 标签）；多个 BAM 的位点合并成一张表。
- velocyto / STARsolo `Velocyto`：没有 velocyto 的两遍 validated intron 逻辑和 exon-intron span 的专门分类（跨边界的块直接算 intronic）；按基因全长判断多基因归属；同一 UMI 的读段落在不同基因时在每个基因各算一个分子（velocyto 丢弃这样的分子）；无 UMI 时 read pair 按 read name 的 64 位哈希归为一个分子。
- 代谢标记（GRAND-SLAM、SLAM-DUNK 一类）：不做读段两端的 trimming、不估计转换率或新旧 RNA 比例（只给计数与直方图），标记判定用固定阈值；分子的标记取其中任一读段，不把多条读段的转换合计；`Labeling` 的直方图每条读段单独计（read pair 的两端各一条）。

- Bismark / MethylDackel：输出单个 BAM 的合并表（不按 OT/OB/CTOT/CTOB 分文件，没有 `bismark2bedGraph` / `coverage2cytosine` 的全基因组胞嘧啶报告，没有覆盖度的胞嘧啶不列出）；重叠 mate 按 Bismark `--no_overlap` 丢 read 2 的重叠部分，而不是 MethylDackel 的 pileup 重叠处理；不提供 M-bias 自动截断（`--ignore` / `--OT` 等参数），M-bias 只报告；没有 `--minDepth` 等深度过滤；Bismark 模式下上下文直接取 `XM`，不复核参考；非 CpG 的 CHG / CHH 不做链合并；`MM` 的 5mC 只看 `C+m`（组合 `C+mh` 中 m 的概率），`h` 等修饰不单独输出，概率阈值固定（不做 modkit 的自动阈值估计）。
- 4DN pairtools：只用主比对，不处理 walk / 嵌合读段（supplementary 排除），没有 pair type 分类（全部 `UU`），没有 `#samheader` 行；去重要求两端位置完全相同（pairtools 默认允许 3 bp 误差）；`#sorted: chr1-chr2-pos1-pos2` 的染色体顺序是 BAM 参考顺序而不是字典序，没有 pairix 索引；自连接过滤按方向和距离而不按限制性片段；`.hic` 不含 expected / 归一化向量，没有 fragment 分辨率。
- cellsnp-lite / vartrix：只计二等位 SNV（多等位、插入缺失跳过并计数）；所有位点都输出为列，不做 `--minCOUNT` / `--minMAF` 过滤；读段单位按 read name 的 64 位哈希去重；默认要求碱基质量 ≥ 20（cellsnp-lite 不按碱基质量过滤）；不做 mode 2 的从头位点发现；UMI 内读段不一致时在每个等位各计一次，不做多数表决；ASE 只用一个样本的定相杂合位点，不做单倍型定相。

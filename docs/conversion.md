# 格式转换

[文档首页](README.md) · [格式与索引](formats.md) · [定量分析](quantification.md)

## 使用

```sh
ow convert genome.fa genome.fa.gz -t 4 -m 1G
ow convert genome.fa.gz genome.2bit -t 4 -m 1G
ow convert cohort.vcf.gz cohort.bcf --index
ow convert reads.sam reads.bam --index
ow convert unsorted.bam sorted.bam --sort --index
ow convert reads.cram reads.bam --reference genome.fa.gz --index
ow convert peaks.bed peaks.bb --chrom-sizes genome.fa.gz.fai
ow convert signal.bedGraph signal.bw --chrom-sizes genome.chrom.sizes
ow convert sample.mcool sample.hic --resolutions 1000,5000
ow convert sample.pairs.gz sample.hic --resolutions 5000,25000
ow convert pileup.bed calls.bb
ow convert pileup.bed 5mC.bw --mod m --value percent --strand +
```

源路径按内容识别，目标路径按扩展名选择；CLI 的 `--from` / `--to` 或 Python 的格式对象可显式指定。`ow convert --help` 列出路由和全部选项。不适用于当前路由的选项会报错。

```python
from omicsweave import convert

plan = convert.plan("cohort.vcf.gz", "cohort.bcf", index=True, threads=4, memory="1G")
print(plan.memory, plan.budget)  # fixed / pipeline / buffer，单位为字节
report = plan.run()  # 或 convert.run(...)：plan 后直接 run
print(report.target, report.companions, report.warnings)
```

`plan` 接受与 `run` 完全相同的选项，检查转换路由、每个选项是否适用和内存规划，不读取数据，也不代替完整输入校验；`Plan.run()` 才处理数据并返回已写出的文件。CLI 的 `--dry-run` 用于查看规划。需要排序的输入必须已排序后才能建索引；SAM → BAM、VCF → BCF 不会隐式排序。

`threads=None` 默认使用进程可用的全部核心（考虑 CPU 亲和性与 cgroup 限制），`memory` 默认 `"4G"`。排序超出缓冲预算时 spill 到输出目录，可用 `tmp_dir` / `--tmp-dir` 指向空间充足的位置。`sidecar=True`（CLI `--sidecar`）同时为 bigBed、bigWig 和 hic 写出 sidecar，默认不写；区间查询不依赖它。BAM `.csi` 的 bin 大小用 `min_shift`，CRAM → BAM 可用 `verify_checksums=False` 跳过校验，bigBed / bigWig 的块与 zoom 布局用 `bbi=BbiOptions(...)`。

所有数据先写入临时文件，成功后发布目标和伴随文件。这不是多个文件之间的文件系统事务；调用方应在转换返回成功后使用整组文件。

BAM → bigWig 的覆盖度、BAM → 计数矩阵属于[定量分析](quantification.md)。

## 转换引擎

`omicsweave-convert` 是端到端的转换工具：在给定的**线程数和内存量**内，流式、并行地把 a 格式文件写成 b 格式文件。

| 转换 | 输入 | 输出 | 对齐基准 |
|---|---|---|---|
| `FastaToBgzipFasta` | `.fa`（纯文本 / gzip / BGZF） | `.fa.gz` + `.fai` + `.gzi` | 与 `compress_fasta`（即 bgzip 的块划分）逐字节一致 |
| `FastaToTwoBit` | `.fa` | `.2bit` | 与 `faToTwoBit` 布局逐字节一致，lib2bit 读回 |
| `VcfToBcf` | `.vcf` | `.bcf`（+ `.csi`） | 未压缩流与 htslib `wbu` 逐字节一致；索引与 `bcftools index` 同一内容 |
| `SamToBam` | `.sam` | `.bam`（+ `.bai` / `.csi`） | 未压缩流与 htslib 逐字节一致；索引与 `samtools index` 同一内容 |
| `CramToBam` | `.cram`（+ 参考） | `.bam`（+ `.bai` / `.csi`） | 未压缩流与 htslib 读 CRAM 写 BAM 逐字节一致（header 为 htslib 解析后重建的文本）；索引与 `samtools index` 同一内容 |
| `BamSort`（`sort`） | `.bam` | 按坐标排序的 `.bam`（+ `.bai` / `.csi`） | 与 `samtools sort --no-PG` 记录顺序、header、块划分一致，level 0 时逐字节一致；索引与 `samtools index` 同一内容 |
| `BedToBigBed` | `.bed`（bedGraph 亦可） | `.bigbed` + `.sidecar` | 与 `BigBedWriter` 逐字节一致 |
| `BedGraphToBigWig` | `.bedgraph` | `.bigwig` + `.sidecar` | 与 `BigWigWriter` 逐字节一致；libBigWig 读回的区间、逐碱基值、总汇总一致，zoom 汇总与区间一致 |
| `WigToBigWig` | `.wig`（fixedStep / variableStep / bedGraph 行） | `.bigwig` + `.sidecar` | 同上 |
| `CoolToHic` | `.cool` / `.mcool` | `.hic` + `.sidecar` | 与 `HicWriter::write_matrix` 逐字节一致 |
| `PairsToHic` | 4DN `.pairs`、Juicer short（纯文本 / gzip / BGZF） | `.hic` + `.sidecar` | 与 `HicBuilder` 逐条 `add` 同样的接触逐字节一致；straw 读回的计数一致 |
| `BedMethylToBigBed` | bedMethyl（modkit / ENCODE） | `.bigbed` + `.sidecar` | 全部字段原样保留，附该方言的 autoSql |
| `BedMethylToBigWig` | bedMethyl（modkit / ENCODE） | `.bigwig` + `.sidecar` | 区间与所选列的值一致；与 `BigWigWriter` 同一写出路径 |

### 选择转换

`Format` 列出 OmicsWeave 读写的全部格式（含 `.fai`、`.gzi`、`.csi`、`.bai`、`.sidecar`），同一格式的不同存储（plain / gzip / bgzip）各是一个值，名称与 Python 的 `File.format` 一致；`Format::detect` 按内容、`Format::from_extension` 按名称识别。`Conversion::between(from, to, src, dst)` 按格式对选择转换，源格式的任一存储都可以；`Conversion::infer` 两端都按扩展名判断；`Conversion::configure(&Options)` 应用选项，不适用于该转换的选项报错。`Report.files` 带每个文件的 `Format`；`Report.warnings` 列出转换跳过而不报错的内容（例如染色体列表外的行），Python 以 `omicsweave.convert.ConvertWarning` 发出，命令行打印到 stderr。

### 分工

Python（`omicsweave.convert`、`ow convert`）只做顶层组装：把路径或文件对象交给原生层，把返回的格式名组装成文件对象，命令行负责参数、进度与输出。识别格式（源按内容、目标按名称）、校验选项、预算缺省值、解析大小都在 Rust。

### 代码结构

`format`、`conversion`（选择、`plan` / `run`、运行上下文）之外，每种源格式一个模块（`fasta`、`vcf`、`sam`、`bamsort`、`cram`、`bed`、`wig`（bedGraph 与 wig）、`bedmethyl`、`cool`、`pairs`（4DN pairs 与 Juicer short）），共用部件在 `pipeline/`：`input`（压缩识别、并行解压、按行分批）、`output`（临时文件与原子改名、spill 目录）、`sort`（外排序）、`records`（VCF → BCF、SAM → BAM 共用的「文本行 → BGZF 记录 + 分箱索引」流水线）、`hic`（cool → hic 与 pairs → hic 共用的「按文件顺序的像素 → 矩阵、块、`.hic` 与 sidecar」）。bedMethyl → bigBed 复用 `bed` 的解析—排序—写出流程（换成 bedMethyl 的行解析与 autoSql），→ bigWig 复用 `wig` 的写出段。

### 三段式流水线

每条转换都由同样的三种阶段组成：

1. **串行读取 / 分帧**：把输入切成单元（行批、BGZF 块、像素区间），并记下只依赖未压缩内容的元数据（`.fai`、zoom 汇总、块号）。开销小。
2. **并行纯计算**：解压、解析、编码、压缩，以任务的形式提交到共享的 `omicsweave_core::exec::Pool`。任务只拥有自己捕获的数据（`'static`），不借用 writer 的状态。
3. **串行按序提交**：按提交顺序取回结果，写出，得到压缩后的偏移，建立依赖偏移的索引（`.gzi`、R-tree、hic 块索引、BAI/CSI、sidecar）。

`core::exec` 的 `Ordered` 把这一接缝标准化：`Inline` 在提交时立即执行，`PoolQueue` 在线程池上执行，结果都按提交顺序返回。各格式的 writer 只保留一条代码路径，`Inline` 只是执行器的一种。由此得到**确定性不变式**：输出与线程数无关，也与排序是否 spill 到磁盘无关，测试对每条转换都检查这一点。

### 哪些能并行，哪些不能


| 转换 | 并行 | 串行 |
|---|---|---|
| FASTA → bgzip | deflate；BGZF 输入的解压 | 读取、`.fai` 扫描、写出 |
| FASTA → 2bit | 各片段打包，统计 N / mask 游程 | 解析 FASTA、合并游程、目录 |
| VCF → BCF、SAM → BAM | 按行批解析编码、deflate；BGZF 输入的解压 | 切行、BGZF 分帧、建索引、写出 |
| BAM → 排序 BAM | 解压、按批校验记录、分区排序、deflate | 切记录、header、归并、BGZF 分帧、建索引 |
| CRAM → BAM | 按 slice 解压各块、解码记录、编码 BAM 记录（参考由共享读取器加锁读取）、deflate | 读容器、定位 slice、BGZF 分帧、建索引 |
| BED → bigBed | 解析、分区排序、块压缩 | 染色体编号、zoom sweep、R-tree |
| bedGraph / wig → bigWig | 解析、分区排序、块压缩 | 解析 wig 段首（lead）、染色体编号、zoom 汇总、R-tree |
| cool → hic | chunk 读取与解码、像素定位、分区排序、块编码与压缩 | 矩阵记录、footer |
| pairs → hic | 解压（BGZF）、按行批解析与最细分辨率分箱、分区排序、按固定大小的像素块汇总到各较粗分辨率、块编码与压缩 | 读 header、切行、合并同一格的计数、矩阵记录、footer |
| bedMethyl → bigBed / bigWig | 解析与校验、按 mod / strand 选择、分区排序、块压缩 | 染色体编号、重叠检查、zoom 汇总、R-tree |

普通 gzip（非 BGZF）无法拆分，只能在单独的线程里流水线解压，不能并行。

### BGZF（`omicsweave-bgzf`）

fasta、vcf、sam、bam 共用。`Writer` 按 htslib `bgzf_write` / `bgzf_flush_try` / `bgzf_flush` 的规则切块，所以未压缩块与 htslib 一致。块在执行器上压缩；写入时报告的位置是 `Logical` 偏移（块序号代替压缩地址），文件写完后由 `Blocks::resolve` 换算成虚拟偏移。`BinningIndex` 按 `hts_idx_push` / `hts_idx_finish` 的逻辑建 BAI / CSI：先在 `Logical` 偏移上做在线合并，偏移换算之后再执行那些需要比较地址大小的步骤（合并小 bin、合并同块 chunk）。`Reader` 在执行器上预读并解压后续块。读取端的 `BinIndex`（载入 BAI / CSI，按 `hts_itr_query` 算出区间的 chunk）、`ChunkWalk`（按 `hts_itr_next` 遍历 chunk）和 `SeekReader`（按虚拟偏移读文本，每个 chunk 一次 pread）由 BAM 与 BCF 的区间查询共用，记录格式各自解析。

### 预算

`Budget { threads, memory, tmp_dir }`。内存采用规划式分配：`plan` 在读数据之前给出 `MemoryPlan`，包括三部分：
- `fixed`：常驻部分，如索引、header、bins 表；
- `pipeline`：每个队列的在途上限 × 单元大小；
- `buffer`：需要重排的转换使用的排序缓冲区。

`fixed + pipeline` 超出预算时直接报错，不在运行中途 OOM。批大小随预算自适应（64 KiB 到 4 MiB）。核算只覆盖大块缓冲，分配器开销、线程栈和小元数据不计入，因此该预算不是进程 RSS 的硬上限，应为分配器、线程栈和运行环境留出余量。

### 超过内存时

- **纯流式**（FASTA → bgzip、VCF → BCF、SAM → BAM）：工作集只有在途批次，任意大小都能处理。
- **单元级**（FASTA → 2bit）：2bit 记录头在打包数据之前，所以一条序列的打包数据先暂存；超出缓冲时写到临时文件。目录在所有记录之前，序列名来自源文件的 `.fai`，没有 `.fai` 时多扫一遍输入。
- **需要重排**（BED / bedMethyl → bigBed、bedGraph / wig / bedMethyl → bigWig、BAM 排序、cool / pairs → hic）：用通用外排序 `pipeline::sort::Sorter`。缓冲满时，并行排序各部分，归并写成 run 文件，最后与剩余的内存部分做 k 路归并。键唯一（BED、BAM 加输入序号，bedGraph / wig 加行号），所以结果与在哪里 spill 无关。BAM 排序至少需要 1 MiB 的排序缓冲区（`samtools sort` 每线程的下限），预算容不下时报错。cool → hic 按 (染色体对, 块号, y, x) 排序，之后用 hic 的流式接口（`begin_matrix` / `add_block` / `end_matrix`）逐块写出，整张矩阵不需要放进内存。pairs → hic 的键相同；同一格的多条接触就是一个像素的计数，`Spill::ABSORBS` 让排序器在写 run 文件时就把相等的键合并（`absorb`），最后的归并流里调用方再合并一次。计数是整数，用 `f64` 求和是精确的，合并发生在哪里都不影响结果；Juicer short 的 score 可以是小数，这时最细分辨率的键带上行号、较粗分辨率的键带上像素块号，排序器不会合并它们，每个格子的和总是按同一顺序求得，与线程数和 spill 无关。

spill 文件默认放在输出文件旁边（不用 `/tmp`）。所有输出先写到 `.partial`，整个转换成功后才改名；失败时清理。

### 内存数据结构

每种格式只定义转换需要的流式单元，不定义整文件模型，也不对外提供操作接口：

| 数据模型 | 单元 |
|---|---|
| VCF / BCF | 复用 `Record` |
| SAM / BAM | 行到 BAM 记录的编码器（`BamEncoder`）；BAM 记录（`Record`） |
| FASTA / 2bit | `TwoBitPiece` / `TwoBitRecord` |
| Interval | `(chrom, start, end, rest)` |
| Contact | 按 `MatrixLayout` 定位后的 `Record` |

查询接口仍然只有"a 文件 → b 内存产物"一种；Python 层只暴露 `omicsweave.convert`（`run` / `plan` / `routes`），命令行是 `ow convert`。

### 与参照工具的有意偏差

- VCF 中 header 未声明的名字报错。htslib 会动态补 header，但写出的 BCF 会损坏；需要补 contig 时Python 开启 `add_contigs=True`（CLI：`--add-contigs`），先扫一遍 CHROM 列。
- SAM 中 RNAME 不在 @SQ 里时报错。htslib 只警告，并把读段当作未比对。
- BED 的染色体按首次出现顺序编号；bigBed 的 B+ 树按名字索引，读取不受影响。bedGraph / wig → bigWig 同样如此（`wigToBigWig` 按名字排序）。
- wig 的 fixedStep / variableStep 段写成 bedGraph 类型的数据块（`wigToBigWig` 保留段类型）；解码出的区间与值相同。bedGraph / wig 不要求按序（UCSC 要求已排序输入），`track` / `browser` 行在任意位置跳过，span 大于 step 造成的重叠报错；空区间、非有限值报错；不提供 `-clip`，越过染色体末端报错。
- wig 按批并行解析：一批开头尚未见到段首的行先记为 lead，由串行阶段按上一批结束时的段补全位置，每行只解码一次。
- cool 只取上三角（`square` 模式同样如此），不读 balance weight。
- `.hic` 不加 Juicer 的 `All` 染色体。
- pairs → hic 不做 Juicer `pre` 的任何过滤（MAPQ、pair type、同一酶切片段内的读段对）：每行算一次接触（short-with-score 算它的 score）。位置按原样分箱（与 `pre` 相同，1-based 的位置不减一），超出染色体长度的位置计在末端（与 `HicBuilder` 相同）。染色体列表外的行跳过，在 `Report.warnings` 里按染色体给出行数。只写 BP 分辨率；FRAG 需要酶切位点文件，不支持。缺省分辨率是 `pre` 的 9 个（2.5M 到 5k）；较粗分辨率由最细的像素汇总，必须是最细分辨率的整数倍，否则报错并指出是哪一个。每格计数先精确求和再存成 `f32`（`pre` 与 `HicBuilder` 直接用 `f32` 累加，计数超过 2^24 后才会不同）。分辨率按 `pre` 的习惯从粗到细列出、按此顺序写出，最细分辨率的像素先写到 spill 目录里的临时文件，最后写出。
- bedMethyl 按字段数分方言：18 列为 modkit `pileup`（第 10–18 列在旧版 modkit 中以空格分隔，同样接受），11 列为 ENCODE WGBS；混用报错。→ bigWig 时值原样保存（百分比保持 0–100，不从计数推导比例）；文件里有多个修饰代码（第 4 列）而没有给 `mod` 时报错并列出代码；所选条目仍有重叠时报错，并指出该加的选项（两条链重叠时提示 `strand`）。
- BAM 排序不加 `@PG` 行（相当于 `samtools sort --no-PG`）；`--index` 缺省写 `.bai`（参考序列放不下时写 `.csi`），`samtools sort --write-index` 缺省写 `.csi`。BAM 到 BAM 只有排序这一种转换，没有 `sort` 时报错；未排序的 BAM 不能只建索引。

### BAM 排序

排序键与 `samtools sort` 的坐标序（`bam_sort.c` 中的 `ks_radixsort` 和归并时的 `heap_lt`）相同：参考序列编号按无符号数比较（`tid` -1 即无参考序列的记录排最后），再比 `pos + 1`（取低 63 位），再比反链标志，其余按输入顺序。记录按原始字节搬运，不重新编码，只写回 htslib 读入时重算的 `bin`（有 CIGAR 的记录）。header 按 htslib 读入、解析、重建的规则改写：`sam_hdr_sanitise`（每行以 `@` 开头，补上最后的换行）、行格式校验、二进制参考序列列表中没有 `@SQ` 行的补在最后一个 `@SQ` 之后、`LN` 不一致时以二进制列表为准、只在文本中出现的 `@SQ` 追加到参考序列列表；第一个 `@HD` 移到最前，`SO` 改为 `coordinate`，删去 `GO`、`SS`，没有 `@HD` 时补 `@HD VN:1.6 SO:coordinate`。htslib 拒绝的 header（缺 `SN` / `LN` 的 `@SQ`、缺 `ID` 的 `@RG` / `@PG`、重复的参考序列等）同样报错。与 `samtools sort` 的差异只在病态输入：read name 缺结尾 NUL、CIGAR 存放在 `CG` 标签中的记录按原样复制（htslib 会补 NUL、把 `CG` 移到标签末尾或把短 CIGAR 写回原位）。参照实现 `refs/samtools`（samtools 1.25 的 `bam_sort.c` 及其依赖，进程内调用）只用于测试。

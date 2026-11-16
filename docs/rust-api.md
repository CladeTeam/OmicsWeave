# Rust 格式接口与存储语义

[文档首页](README.md) · [Python 查询](queries.md) · [架构](design.md)

本文记录 Rust 层的访问契约、布局和格式兼容边界。Python 用户从[查询指南](queries.md)开始；Rust 符号的完整签名可在仓库运行 `cargo doc --workspace --no-deps` 后查看。`refs/` 仅用于差分测试，不是公共接口。

## 标准访问与输出

同一系列的不同格式，访问方式和输出布局完全相同；格式之间的差异只存在于各自的低层步骤（算字节范围、解码块）中，这些步骤保留给自行调度 I/O 的 loader 使用。

### 访问 key

坐标统一为 `u32`，0-based、左闭右开，与 BBI、2bit、bigBed 等格式的存储宽度一致；文件偏移等元数据仍为 `u64`。key 分两种：

- `Region { chrom, start, end }`：一条染色体上的一个区间，是稀疏区域查询（特征、读段、记录）的 key；contact 使用两个 Region 组成的 ContactQuery。
- `Locus { chrom, spans, strand }`：同一条染色体上升序、互不重叠（可以相接）的若干片段，加一条链 `Strand`（`+`、`-`、`.`、`?`），是稠密查询（碱基、逐位信号、bin 统计）的 key。`Region` 可直接转成单片段、无链的 locus；剪接后的转录本就是由外显子组成的 locus。

### 方向与拼接的统一规则


- **稠密输出**按 locus 的 5'→3' 写出：各片段首尾拼接，`-` 链整体翻转，碱基再取互补（IUPAC 码映射到互补集合的码，保留大小写，`U` 按 `T` 处理）；`.` 和 `?` 按 `+` 读。长度 = 各片段长度之和。
- **稀疏输出**带自己的坐标：保持基因组坐标和文件顺序，strand 既不过滤，也不翻转。
- **只解码一次**：先求出与任一片段重叠的数据块并去重，每块只读、只解压一次，再分发到拼接后的输出位置。不会去解码片段之间的内含子；多个片段共享同一个块时也不重复解码。bgzip FASTA 中，共享 BGZF 块的相邻片段合并成一段一起解码。
- **翻转在原地进行**：先按正向解码进调用方的 buffer，再原地 reverse（`core::orient`），不需要额外内存。

### 五个系列

（另有 BCF 的基因型矩阵，见下文「变异（BCF）」）：

| 系列 | 格式 | 访问 trait | 输出 |
|---|---|---|---|
| Sequence | `.fa`、`.fa.gz`、`.2bit` | `SequenceRead` | `[u8]` ASCII，长度 = `locus.len()` |
| Track | bigWig | `TrackRead` | values：`[f32]`，长度 = `locus.len()`，无数据为 NaN；bins：`[f32]`，每个统计量一行（见下文「bin 汇总」）；entries：`[TrackEntry]`，**裁剪到查询范围** |
| Interval | bigBed | `IntervalRead` | `[IntervalEntry]` + 注释字符串列；**不裁剪** |
| Contact | `.hic` | `ContactRead` | 稠密 tile：行优先 `[f32]`，形状 = 行 bin 数 × 列 bin 数；或 `[ContactEntry]` 非零格（按 (row, col) 排序） |
| Alignment | BAM（+ `.bai` / `.csi`）、CRAM（+ `.crai`） | `AlignmentRead` | `[AlignmentEntry]` + 按 `AlignmentFields` 选取的变长列；**不裁剪**，文件顺序 |

- **行优先的定长 entry**：`TrackEntry { start: u32, end: u32, value: f32 }`、`IntervalEntry { start: u32, end: u32, annotation: u32 }`，都是 12 字节的 `#[repr(C)]` 结构，与 numpy 结构化 dtype 一一对应。
- **注释字符串列**：Arrow `LargeString` 布局（`i64` offsets 长度 n+1，加 `u8` data）。`IntervalEntry::annotation` 是该列中的行号，所以对 entry 排序时不必移动字符串。每行内容是记录第三列之后的全部字段，用 tab 连接，与 bigBed 的存储方式相同。
- **Track 裁剪、Interval 不裁剪**：Track 的区间只是信号的一种表达方式，不应影响数据内容，因此查多少给多少；Interval 的每条记录是完整的特征（例如一个 peak），注释描述的是整个区间，因此保持原样。
- **变长输出：一次解码，暂存，再拷贝**：`TrackRead::query_entries` / `IntervalRead::query` 对查询区间只解码一次，把结果暂存在 reader 内部并返回输出大小；调用方分配好内存后，`write_entries` / `write`（`&self`，不再读取或解码）把暂存结果拷贝进调用方的缓冲区。暂存属于 reader 的内部工作状态，交付的输出始终在调用方的内存中。缓冲区大小不符返回 `Error::OutputSize`，未查询就写返回 `Error::NotQueried`。
  - bigBed：查询时把解压后的块依次追加到一块只增不缩的内部缓冲区，暂存的记录只保存注释在其中的位置，不拷贝注释字节；`write` 是唯一一次拷贝。
  - bigWig entries：暂存解码并裁剪好的 `TrackEntry`，`write_entries` 是一次 memcpy。
- **越界与未知染色体**：Sequence 的输出必须定长，任一片段超出序列即报 `RegionOutOfRange`；Track 和 Interval 允许超出染色体末端，没有数据的部分就是没有。未知染色体：Sequence 和 Track 报 `UnknownChrom`；Interval 返回空。空区间在所有系列中都返回空。
- **统一 Error**：所有格式 crate 使用 `omicsweave_core::Error`，Python 绑定映射同一套异常。
- **字节源**：读取器通过调用方提供的 `ReadAt`（`[u8]`、`Vec<u8>`、`FileSource` 等）读取所需字节；I/O 来自何处、是否缓存，由调用方决定。

### 不重复解码

任何重复扫描或重复解码，只有在开销小于一次拷贝时才允许；否则一律解码一次、暂存、再拷贝。定长输出（Sequence、Track values）直接解码进调用方的缓冲区，不经过暂存。

### bin 汇总（bigWig）

`TrackRead::bins(locus, BinSpec { resolution, phase, stats, source }, out)` 在读取器内把信号汇总成定宽 bin，不展开逐碱基数组。

- **bin 的划分**：bin 沿 locus 的 5'→3' 排列，局部位置 `x` 落在 bin `(x + phase) / resolution`，共 `ceil((len + phase) / resolution)` 个。`phase = 0` 且长度是 `resolution` 的整数倍时每个 bin 都完整；否则首尾 bin 不完整。如何对齐由调用方决定：例如正链单片段 locus 要对齐到基因组网格，就传 `phase = start % resolution`。反链在内部换算成镜像的正向 phase，算完后把每行翻转。
- **统计量**：mean、min、max、sum、std、coverage，可以一次要多个，输出 `[统计量数, bin 数]`，按请求顺序排列。所有统计量共用一次解码、一次扫描；std 的第二遍只扫已暂存的数据段，不重新解码。
- **语义**：与 libBigWig 的 `bwStatsFromFull` 一致，完整 bin 上逐位相同：
  - mean 是对有数据碱基的长度加权平均；
  - sum 与 libBigWig 相同，先用 `float` 计算 `len × value` 再累加；
  - std 是两遍算法的样本标准差，只有一个碱基时为 0；
  - min/max 取与 bin 重叠的数据段；
  - coverage = 有数据的碱基数 / bin 内的碱基数；
  - 全部用 `f64` 累加，最后舍入成 `f32`。
- **空 bin**：除 coverage 外都是 NaN；coverage 为 0。libBigWig 在这里给 NaN，这是一处有意的差异。「缺失按 0 算的均值」等于 sum / bin 宽度。
- **zoom**（`BinSource::Zoom`）：按 libBigWig `bwStats` 的规则选级别（缩减倍数不超过半个 bin 的最粗级别），每条 zoom 记录按落在 bin 内的比例计入；没有合适级别时回退到精确计算。mean、min、max、std、coverage 与 `bwStats` 逐位一致。sum 是有意的差异：libBigWig 把比例截断成整数，我们按比例加权。zoom 只是估计值，默认用精确计算。
- **对齐基准**：差分测试随机生成 bedGraph / varStep / fixedStep 混合的 bigWig，与 `bwStatsFromFull` / `bwStats` 逐位比较；对任意多片段、带链、带 phase、不完整 bin 的 locus，与逐碱基 values 的暴力汇总比较。

### Contact（Hi-C）

只支持 `.hic` v8 和 v9，其它版本报 `Error::Unsupported`。只提供原始接触计数：归一化和 expected 属于分析逻辑，格式层不读取，也不计算。查询是 `ContactQuery { rows, cols, resolution, unit }`，行列均为连续 `Region`，不接受多片段 `Locus`。Python 的 `ContactWindow` 另行处理两条轴的 strand 与输出翻转。

- **区间与方向**：行和列可以在不同染色体上，区间覆盖的每个 bin 都计入 tile。无论文件存的是对称矩阵的哪一半，都按查询方向给出完整的 tile（对角线另一侧自动镜像）。
- **单位 `unit`**：`Bp` 时区间是碱基坐标；`Frag` 时区间是酶切片段编号，片段 `i` 从第 `i-1` 个酶切位点到第 `i` 个。
- **稠密输出**：没有接触的格为 0。
- **entries**：只含计数有限且非零的格，坐标为 bin 起点（单位同查询）。
- **只解码一次**：块索引和酶切位点首次使用时读入并缓存；一次查询中每个块只解压、解码一次。
- **v8 与 v9 的差别**只在四处：header 中的 NVI 指针和染色体长度宽度、footer 中 `nBytesV5` 的宽度、块内坐标宽度标志（v9 才有）、染色体内的块编号（v9 沿对角线按深度编号，受 `v9-depth-base` 控制）。
- 读出的计数以 straw 的 `observed` / `NONE` 为参照逐位校验。

### 写 `.hic`

- `HicWriter` 按 Juicer 的布局写 v8 或 v9（默认 v9）：header、blocks、矩阵记录、footer、normalized expected、norm vector index。expected 和归一化列表写为空，straw 仍可用 `NONE` 读取。
  - 调用方给出分好 bin 的矩阵。
  - 每个 block 选用最小的表示：list-of-rows 或 dense，计数用 short 或 float；v9 中坐标可在 short / int 之间按需选择。
  - 支持 v9 的 `v9-depth-base` 布局。
  - **块大小可配置**（`BlockLayout` / `BlockBins`，`HicWriter::write_matrix_with_blocks` 可按单个矩阵覆盖）。默认染色体内矩阵用 128 bin 的块，染色体间矩阵沿用 Juicer 网格（约 1000 bin）。理由：按接触密度 ∝ 1/距离 估算，1 kb 下一个 131 kb 的 cis 窗口在 v9 对角线布局中，用 1000 bin 块要解码约 15 倍于所用的接触，用 128 bin 块约 2.5 倍、只涉及约 2 个块；更小的块主要增加请求数。染色体间矩阵稀疏，块切小只会产生大量碎块。
  - 块大小写在每个矩阵记录里，任何读取器都按文件中的值解析；straw 对 1–400 bin 的块与我们逐位一致。注意 straw 计算 v9 染色体内块号时固定用以 2 为底的 log 深度，忽略 `v9-depth-base`，所以非默认深度布局的文件只有我们的读取器能正确读。
- `HicBuilder` 从接触对构建完整文件，相当于 Juicer 的 `pre`：生成所有 BP / FRAG 分辨率，不计算 expected 和归一化。
- 写出的文件以 straw 为参照逐位校验。

### 个体化序列（consensus）

`omicsweave-consensus` 是 `bcftools consensus` 1.23 的移植，源自 pyconsensus。
- **输入与输出**：只读 BGZF BCF 与它的 CSI，经 `omicsweave-vcf` 的 `BcfReader` 按区间读取（见下文「变异（BCF）」）；VCF 文本先 `ow convert x.vcf.gz x.bcf --index`。`load` 只解码索引给出、且与区间重叠的记录，每条一次；每条记录只留构建要用的部分：位置、参考跨度、等位，以及按 bcftools 读法保存的 `int8` 基因型。取哪些记录与 htslib synced reader 对 `--regions-overlap` 的判断一致：先按 `bcf_itr_queryi` 的存储跨度取候选，再按解码出的跨度（`0` 看起点，`1` 看整条记录，`2` 看去掉与所有 ALT 共同前缀后的部分）筛选，所以起点在区间之前的删除等记录也会被取到。参考序列经任意 `SequenceRead` 读取（`Locus::from(region)`）。输出长度随 indel 变化，沿用"查询返回大小、再拷进调用方缓冲区"的模式。
- **链**：`write(strand, out)` 按链给出序列：`-` 时边拷贝边反向互补（IUPAC 码取互补集合，保留大小写），其余链按正向输出；chain 始终是正链坐标。Python 层按 `Region.strand` 传入，批量结果中每条序列在自己的槽位内反向；`Spliced` key 不适用（consensus 只编辑一段连续区间），报 `TypeError`。
- **批量**：`load` 读一次区间的参考序列（含 mask），并解码该区间的记录；之后每次 `query` 按 `Alleles` 生成一条序列（对应 `-s`、`-H`、`-I`）。同一区间的大量样本、单倍型共用这一次加载，每条序列的代价只剩线性拷贝和逐记录的等位选择。
- **单一执行路径**：逐条移植 `apply_variant` 的状态机，不设快速路径分支。consensus.c 原地修改缓冲区，每遇到 indel 就整体平移尾部；这里改为"已改写的前缀 + 未改动的参考后缀"，记录按序到达，且不会改写前一条记录之前的位置，所以每个区间的耗时与长度成线性。
- **支持的选项**：`bcftools consensus` 的选项全部支持，包括 `-I`、`-M`、`-a`、`-m`/`--mask-with`、`--mark-*`、`-c`、`--regions-overlap` 和 `-i`/`-e`。
- **过滤**：consensus.c 对不满足 `-i`/`-e` 的记录只是跳过，所以过滤在 `load` 解码记录时做（`omicsweave-vcf` 的 `Filter`，即 `filter.c` 的重写），按文件顺序逐条求值；表达式在构建 `Consensus` 时编译一次，每个区间从未求值过的状态开始，与只含该区间的一次 `bcftools consensus` 一致，结果也不依赖查询顺序。表达式编译失败时构建报错；某条记录求值出错，则记在该记录上，查询经过它时才报，与 bcftools 读到该记录才报错的时机一致。
- **不发明缓存格式**：不使用 pyconsensus 的 `.cvcf`，也不把整个文件读成常驻对象：BCF + CSI 已能按区间只读所需的块。克隆（`clone_context` / `try_clone`）共享 header、索引、编译好的过滤式和选项，各自持有读取与解码状态。
- **对齐基准**：`refs/bcftools` 在进程内运行原版 bcftools，以它的输出为准逐字节对齐，包括它的若干行为怪癖。

### 比对（BAM）

`omicsweave-bam`，以 htslib 的读法为准。
- **打开**：读 header，记录的参考编号以二进制参考列表为准；参考名必须是 UTF-8，重名时按名字找到第一个。索引（BAI 或 CSI）由调用方给出；Python 层按 htslib 的顺序在旁边查找（`<bam>.csi`、`<name>.csi`、`<bam>.bai`、`<name>.bai`），找不到即报 `IoError`，打开时从不构建。
- **查询**：按 `hts_itr_query` 由索引算出 chunk（linear index 与 `loff` 给出 `min_off`，右侧 bin 给出 `max_off`，再去包含、去重叠、合并同块），按 `hts_itr_next` 的规则读：每个 chunk 一次 pread，记录的 `[pos, bam_endpos)` 与区间重叠才返回，遇到其它参考或 `pos ≥ end` 即停。结果与 `sam_itr_queryi` 逐条、逐字段一致，包括放在参考上的 unmapped 记录。随机读不依赖 BGZF 之外的派生产物。
- **记录**：校验与 `bam_read1` / `sam_read1` 相同（长度、参考编号范围、CIGAR 与序列长度）；超过 65535 个操作的 CIGAR 放在 `CG:B:I` 标签里时按 `bam_tag2cigar` 还原，aux 输出不含该标签。
- **输出**：`AlignmentEntry`（32 字节：`start` / `end` 为 `i64`，`end` 即 `bam_endpos`；`mate_start`、`mate_tid`、`tlen` 为 `i32`；`flag`、`mapq`），加上按 `AlignmentFields` 选取的变长列（name、ASCII 序列、quality、`u32` CIGAR、原始编码的 aux），都是 offsets + data。查询时只把选中字段的字节暂存一次，`write` 解包碱基并拷出。
- **建索引**：`build_index` 顺序读整个文件（解压并行），按 `bgzf_tell` 的语义取每条记录的结束位置（记录正好在块尾结束时，位置指向下一个物理块，空块也算），交给与 `SamToBam` 共用的 `BinningIndex`，结果与 `samtools index` 内容相同。
- **对齐基准**：htslib 的 `sam_itr_queryi` 与 `sam_index_build3`（`refs/htslib/src/bam.c` 导出 `bam1_t` 字段和 `sam.h` 中的内联迭代器），覆盖 BAI、不同 `min_shift` 的 CSI、多种压缩级别、文件中间的空 BGZF 块和长 CIGAR。
- **不做**：SAM 不支持随机访问，只作为转换源；不做覆盖度、pileup 或从读段生成 consensus：计数、变异检测属于业务逻辑。

### 比对（CRAM）

`omicsweave-cram`，以 htslib 1.25 的读法为准，经 `AlignmentRead` 给出与 BAM 完全相同的输出。
- **版本与解码**：只读 CRAM 3.0 与 3.1，其它版本报 `Error::Unsupported`。容器、slice、压缩头和各数据序列的编码（EXTERNAL、HUFFMAN、BETA、SUBEXP、GAMMA、BYTE_ARRAY_LEN、BYTE_ARRAY_STOP）用 Rust 解码；块的 gzip 用 libdeflate，bzip2、lzma 用静态链接的 `bzip2-sys`、`lzma-sys`，rANS 4x8 / Nx16、自适应算术编码、fqzcomp 与名字分词器来自 `htscodecs-sys`（htscodecs 1.6.7 原样 vendor，SIMD 解码器按 CPU 运行时选择）。块的 CRC32 在解压时校验，从不解压的块不校验，与 htslib 相同。
- **参考序列**：调用方给出任意 `SequenceRead`（bgzip FASTA、2bit 等），按 `@SQ` 名字取 slice 的参考区间并转大写；与 htslib 一样校验 slice 的参考 MD5（可与 CRC32 一并关闭）。内嵌参考与无参考（`RR=0`）的文件不需要参考；多参考 slice 只读取查询所在参考的序列。
- **索引**：`.crai`（gzip 文本）按 `cram_index_load` 载入成嵌套包含表，起点由 `cram_index_query` 的二分查找给出；Python 层在旁边找 `<cram>.crai`、`<name>.crai`，缺失即报 `IoError`，打开时从不构建。`build_index` 只解码多参考 slice（按记录求各参考的范围，同 `cram_index_build_multiref`），其余条目直接取自 slice 头，解压后与 `samtools index` 逐字节一致。
- **查询**：从索引给出的容器顺序读，按 `cram_first_slice` / `cram_next_slice` 跳过或终止容器与 slice，再按 `cram_get_seq` 用记录的 1-based `apos` / `aend` 过滤；结果与 `sam_itr_queryi` 逐条、逐字段一致，包括多参考 slice、放在参考上的 unmapped、超长 CIGAR 与 `CG` 标签。空区间返回空：htslib 对 CRAM 会返回覆盖该位置的记录，这里与 BAM 和 `AlignmentRead` 的约定一致。
- **记录**：按 `cram_decode_slice` / `cram_decode_seq` / `cram_to_bam` 还原：slice 内 mate 的交叉引用与 TLEN 推算，没有存名字的读段按文件名和记录号命名，读组追加为 `RG:Z`，未存储的 MD / NM 按参考生成（htslib 默认如此），`CG:B:I` 按 `bam_tag2cigar` 移回 CIGAR。碱基经 4-bit 编码往返，与 BAM 读出的一致。
- **按字段解码**：CRAM 的各数据序列分块存放。查询由所需字段推出数据序列集合，再闭包到与之共享块（或 CORE 比特流）的序列与标签，只解压用到的块。只要定长字段与 CIGAR 时，解码 BF、CF、RI、RL、AP、MQ、mate 相关序列与读段特征的位置和长度，由特征推出 CIGAR 和比对终点，不重建碱基，也不需要参考；碱基、质量值、名字与标签所在的块都不解压。htslib 完整解码时，持有参考的 slice 中越过参考末端的读段终点算法不同；这里按"完整解码是否持有参考"保留这一差别，所以任意字段组合的结果都与完整解码一致。名字、质量值同样不需要参考；碱基与 MD / NM 需要。
- **对齐基准**：`refs/htslib`（启用 bzip2、lzma）在测试中按 `hts_opt_add` 的格式选项写出 3.0 / 3.1、九种块压缩方法、内嵌参考、无参考、多参考 slice、小 slice 和多 slice 容器的 CRAM，与 `sam_itr_queryi`、`sam_index_build3` 逐字段、逐字节比较。

### 变异（BCF）

`omicsweave-vcf` 的 `BcfReader`，以 htslib 的读法为准。

- **打开**：只读 BGZF 压缩的 BCF 与它的 CSI；VCF 文本不按区间读，先转换成 BCF（`ow convert --index` 同时写出 `.csi`）。索引由调用方给出，Python 层在旁边找 `<bcf>.csi`，缺失即报 `IoError`。
- **查询**：与 BAM 共用 `omicsweave-bgzf` 的 `BinIndex`（BAI / CSI 的载入与 `hts_itr_query` 的 chunk 计算）、`ChunkWalk`（`hts_itr_next` 的 chunk 遍历）和 `SeekReader`（按虚拟偏移读 BGZF 文本）；记录的 `[pos, pos + rlen)` 用存储的 `rlen` 判断重叠，与 `bcf_itr_queryi` 逐条一致。只解码重叠的记录。
- **基因型矩阵**：`genotypes(region, samples)` 一次给出区间内每条记录的位置、参考跨度（htslib 解码的 `rlen`）、ID、REF、逗号连接的 ALT（Arrow 字符串列）和所选样本的 `int8 [记录, 样本, 倍性]` 等位基因编号（`-1` 缺失、`-2` 无此等位基因），以及 `bool [记录, 样本]` 的 phase。倍性取区间内最大值；各记录先按自身倍性暂存，倍性都相同时直接作为结果，否则只重排一次。剂量、频率等由调用方从矩阵计算。
- **对齐基准**：htslib 的 `bcf_itr_queryi` / `bcf_itr_next` 与 `bcf_get_format_values`，覆盖 SNP、indel、多等位、带 `END` 的符号等位、混合倍性、缺失和无 GT 的记录，以及不同 `min_shift` 的 CSI。

## Sidecar

loader 要在读取数据块之前静态分配内存，因此需要知道每个中间产物的大小。有些格式没有记录块解压后的大小：bigWig / bigBed 的 R-tree 和 `.hic` 的块索引只给出压缩字节数，zlib 流本身也不记录原始长度；`.hic` 的块编码（list-of-rows 或 dense，short 或 int）还是按块选择的。sidecar 为每个块记录这些信息，属于可重建的派生产物：可以随时删除，并按源文件身份重建。

- **格式**：紧凑的定长小端二进制，不压缩。`omicsweave-core::sidecar` 定义公共的 48 字节容器头：magic `OWSIDECR`、容器版本、格式种类、body 版本、源文件身份、body 长度和校验和。body 由各格式定义，由 8 字节对齐的定长记录表组成，读入后可以直接按下标访问。
- **源文件身份**：文件长度，加上对开头 64 KiB、结尾 64 KiB 和 14 个均匀分布的 4 KiB 页做的指纹。它只依赖文件内容，不依赖路径或 mtime，所以源文件和 sidecar 一起复制、移动后仍然有效。身份不符时报 `Error::StaleSidecar`，调用方应据此重建。
- **bbi**（`BbiSidecar`）：base 级别和每个 zoom 级别的所有块，每块记录 span、文件偏移、压缩大小、解压大小和条目数。块在级别内的位置就是它的稠密序号。`BigWigWriter` / `BigBedWriter` 在写文件时同时产出 sidecar，`BbiSidecar::build` 从已有文件构建。bigWig 也可以不用 sidecar，按头部的 `uncompressBufSize` 上界分配。
- **hic**（`HicSidecar`）：每个 (chr1, chr2, unit, bin size) 矩阵的所有块，每块记录块号、位置、压缩大小、解压大小和解码出的记录数。`HicWriter` / `HicBuilder` 的 `finish_with_sidecar` 在写文件时同时产出 sidecar，`HicSidecar::build` / `build_filtered` 从已有文件构建（后者只处理选中的矩阵）。

## FASTA 索引

`.fa` / `.fa.gz` 的索引沿用标准格式，与 `samtools faidx` / `bgzip` 逐字节一致，接口与 sidecar 同一风格：`Fai` / `Gzi` 各有 `build`（从源文件构建）、`parse` / `to_bytes`、`save` / `load`（`<源文件>.fai`、`<源文件>.gzi`，原子写入）和 `build_file`。

- `build_index` 对 bgzip 文件只解压一遍，同时得到 `.fai` 和 `.gzi`；`Gzi::build` 只读每个块的头部和 ISIZE，不解压。
- `FastaWriter` 写文件时直接按写出的布局记下 `.fai` 条目（压缩时同时记下 `.gzi`），不重新解析写出的内容。
- 空的 BGZF 块（包括 EOF 块）不进 `.gzi`，与 htslib 一致。

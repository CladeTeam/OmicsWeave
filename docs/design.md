# OmicsWeave 设计

OmicsWeave 是一个系统性的 a2b 框架：a 端是组学数据的输入格式，b 端是模型训练使用的数据形式。框架负责把任意受支持的 a 端数据组织成 b 端训练数据，并把多种组学数据整合到统一的加载器中。

## 1. 设计原则

### 1.1 沿用现有格式，只在真实瓶颈处设计新格式

a 端格式按原样读取，不要求用户先把数据转换成 OmicsWeave 自有的存储格式。很多生信格式效率不高，例如普通 gzip 压缩的 FASTA 无法随机访问，压缩收益也有限，但只是格式不理想不构成重新设计的理由。

只有当某个 a2b 路径在真实负载下出现可测量的瓶颈，且索引、缓存等手段无法解决时，才设计新的存储格式。CellPack 就是这样产生的：单细胞表达矩阵的 shuffle 训练加载被 H5AD 的读取方式卡住，才有了逐细胞可独立解码的新格式。

由此得到以下约束：

- **原始文件始终是事实来源。** 框架不改写、不替换用户的原始数据。
- **允许派生产物，但它们只是加速手段。** 索引（`.fai`、`.tbi`、`.csi`）、解压副本、编译缓存（例如 pyconsensus 的 `.cvcf`）都属于派生产物：可以随时删除、按需重建，并按源文件身份（路径、大小、mtime 或内容哈希）校验是否失效。派生产物不得成为唯一的数据来源。
- **低效格式的处理顺序**：先用标准索引（如 bgzip + fai）；其次构建派生缓存（如把 `.fa.gz` 解压一次成可随机读取的副本，或常驻内存）；仍不够时，才进入新格式设计。
- **新格式是独立项目。** 像 CellPack 这样的新格式应该能单独发布、单独使用，OmicsWeave 通过依赖接入，不把格式定义藏在框架内部。

### 1.2 进程内链接工具库，不调用子进程 CLI

需要借助现有生信工具（htslib、bcftools、samtools、bedtools 等）的能力时，优先以库的形式静态链接进进程（C 库通过 FFI，或使用成熟的 Rust/Python 原生实现），一般不通过子进程调用其命令行。

子进程方式会把框架重新拉回传统工具链的约束：输入输出必须经过文件系统，每次调用都要付出进程启动和重复解析的开销，数据无法按需在线生成。pyconsensus 解决的正是这个问题：`bcftools consensus` 的逐调用子进程模式无法在线构造个体基因组序列，移植到进程内并常驻 ref 和 VCF 后，才能在训练时按需生成。

由此得到以下约束：

- **第三方 C 源码放在 `vendor/`**，从源码构建并以 `-fPIC` 静态链接进扩展模块，wheel 不依赖系统动态库或 PATH 中的命令。
- **需要修改或移植的工具放在 `forks/`**，保留上游源码作为行为对照；移植实现需要与上游 CLI 做逐字节对拍测试（参照 pyconsensus 与 bcftools 1.23 的 parity 测试）。
- **只把 CLI 当作参照实现和测试基准**，不把它当作运行时依赖。
- **引入工具前检查许可证**，确认其可以静态链接并随 wheel 分发。

### 1.3 在线优先，物化可选

b 端默认按需在线产出：训练时从原始文件（加上派生加速产物）直接构造样本，不要求预先把整个训练集展开写盘。把 b 端数据物化成文件（例如分片、CellPack 文件）是一种可选的输出，用于在线构造成本确实过高的场景。

### 1.4 不使用 mmap，I/O 与预取显式管理

mmap 在生信领域看起来方便，可以直接把 FASTA 或矩阵当成数组访问，实际却很低效：

- 缺页是同步阻塞的：每次随机访问都可能让一个线程停下等 I/O，读请求无法批量提交，也无法与计算重叠。
- 内核 readahead 按局部顺序猜测访问模式，与训练时的 shuffle 访问不匹配，结果要么多读，要么没读到需要的页。
- 淘汰由内核 LRU 决定，内核不知道未来的访问，因此用不上「访问序列已知」这一关键信息。
- 页缓存计入容器 cgroup 的内存，多进程 DataLoader 各自映射同一文件，内存占用和性能都难以预测。
- 在共享或网络文件系统上，每次缺页都可能变成一次网络往返，上述代价被放大。

由此得到以下约束：

- **所有数据读取使用显式 I/O**（io_uring / preadv，需要时用 `O_DIRECT`），读进框架自己管理的缓冲区（frame）。
- **预取由框架根据已知访问序列主动发起**；缓存哪些数据、何时淘汰由 planner 决定，不依赖 OS 页缓存的启发式。
- 派生产物（例如解压后的 FASTA 副本）同样通过显式 I/O 读取，不用 mmap 访问。

### 1.5 组件不持有数值数据，在外部分配的 buffer 上操作

所有工具包（crate）都不负责数值数据的生命周期：输入和输出都是调用方分配好的 buffer，组件只在上面读写。组件不把「自己分配并返回数据」作为一等接口，也不用这种接口替代外部分配。

代价是调用方要多做几步：先读元数据，按元数据算出需要的内存并分配，再自己完成 I/O。换来的好处是：

- **组合性**：组件不关心数据从哪里来、放在哪里，同一个解码函数既能写进运行时的 cache frame，也能写进输出 ring 的 slot，或者测试里的一个普通数组。
- **跨语言简单**：Python 绑定可以直接把 numpy 数组的内存交给组件，结果原地写进 numpy，不涉及跨语言的所有权和生命周期问题。
- **内存可预算**：内存全部由调用方分配，运行时可以精确统计和限制。

这条约束针对数值数据（信号值、序列、解压后的数据块、索引字节等）。轻量的元数据（文件头字段、染色体名表、zoom 级别表等）允许由组件自己持有。

由此得到以下约束：

- **组件接口分成「描述」和「操作」两步。** 描述步骤由元数据算出需要读取的字节范围和需要分配的大小；操作步骤只接受调用方提供的 buffer。
- **解码结果是借用调用方 buffer 的视图，或者直接写进调用方的输出 buffer**，不返回组件自己分配的容器。
- **可复用的工作状态（例如解压器）由调用方创建并持有**，组件不在内部缓存。

## 2. Loader 架构：N + 1 + M

### 2.1 适用场景

CellPack 的 loader 模块虽然随 cpk 格式一起开发，但它解决的问题并不依赖 cpk。它适用于同时满足以下条件的场景：

- **可随机访问**：每个数据单元（item）都能定位到它依赖的字节范围或其他对象；
- **访问序列已知**：要读哪些 item、按什么顺序读，在 run 开始前就确定，或者可以由确定性函数按区间计算（例如按 seed 生成的 shuffle）；
- **可缓存**：重复访问或相邻 item 会共享底层数据，缓存能省下 I/O 或解码开销。

训练数据加载正好满足这三点。由于访问序列已知，I/O 和缓存可以离线或在线（带 lookahead）提前编排，做到 mmap 加 OS 启发式无法做到的预取和淘汰（见 1.4）。

如果每种 a 端格式都配一个面向每种 b 端训练格式的加载器，就是 N × M 的组合。OmicsWeave 把加载拆成三部分，得到 N + 1 + M：

| 部分 | 数量 | 职责 | 不负责 |
|---|---|---|---|
| **Source 适配器** | N（每种 a 端格式一个） | catalog：把文件组织进全局 item id 空间，给出 item 依赖的对象（文件字节范围、解码块等）及其尺寸；decode kernel：把读入的字节解码成该 source 的规范行 | 调度、缓存、I/O、输出 |
| **通用运行时** | 1 | Accesses（position → item 的纯函数）；Planner（根据已知序列、对象尺寸和代价模型决定缓存驻留）；Executor（worker 池，显式预取、等待 frame、解码、映射）；Output ring（固定 slot、多 schema region、lease 背压、多 rank 分发） | 任何格式细节、训练语义 |
| **Output schema / process** | M（每种 b 端训练格式一个） | 把规范行转换成固定大小的输出记录：dense 向量、one-hot 序列、token 序列等 | 数据从哪里来 |

### 2.2 CellPack loader 的现状

| 模块 | 是否通用 | 说明 |
|---|---|---|
| `loader/accesses` | 通用 | 只处理 id 序列：Shuffle、Replicate、Offset、外部回调 |
| `loader/ring`、`planring`、`exec/progress`、`memory` | 通用 | 输出 ring、计划 ring、进度窗口、大页缓冲 |
| `loader/exec/io` | 通用 | io_uring / preadv 后端 |
| `loader/process` | 基本通用 | C ABI 的 process 回调；输入目前是 CSR / dense 行 |
| `opt`（planner） | 算法通用，数据模型耦合 | 复用区间、frame 分配、在线 lookahead 是通用的；`opt::records::Record` 带有 chunk、编码方式等 cpk 字段 |
| `loader/catalog` | 耦合 cpk | `DataSource` 只能注册 `.cpk`，item → (fd, offset, pages) |
| `loader/exec/worker` 的解码步骤 | 耦合 cpk | 直接调用 `scans::PagedRecord` |
| `loader/convert` | 耦合矩阵语义 | 列映射和 dtype 转换，适用于矩阵类 source |

推广的工作主要是把这三处耦合点抽象成 Source 接口，并把 planner 的数据模型从 cpk 记录泛化成通用对象。

### 2.3 推广时要解决的问题

1. **多级派生对象。** CellPack 有两级缓存对象：原始页和解码行。组学格式需要更多层级，例如 BGZF 的原始页 → 解压块 → 解析出的记录 → item 输入。Planner 的对象需要泛化成「带 key、尺寸和生成代价的对象」组成的依赖链。
2. **item 依赖多个对象，甚至跨多个 source。** FASTA 窗口跨越多个页，相邻窗口互相重叠；一个 consensus item 同时依赖 ref 区间、VCF 区间和常驻的编译结构。
3. **常驻对象。** 有些 source 的索引或编译结构整体常驻内存（fai、VCF 编译结构、基因注释），不参与 frame 调度，但要计入内存预算。
4. **变长输出。** 输出 ring 的每个 schema 使用固定的 `cell_bytes`。变长序列、CSR 行等需要 padding 加长度，或者设计 ragged region。
5. **多组学对齐（Weave）就是 Accesses。** 一个 position 的 S 个 access 可以分别指向不同 source 的 item，CellPack 的多 schema 已经支持每个 schema 取不同的 cell。因此 sample / cell / gene / region 之间的 join 归结为构造访问函数，不需要单独的执行路径。
6. **不可随机访问的格式。** 普通 gzip 无法给出 item 的字节范围，需要派生索引（gzip checkpoint 索引）或整体常驻的解压结果，符合 1.1 的处理顺序。

CellPack 与 pyconsensus 是已经成熟的独立项目，OmicsWeave 只借鉴它们的设计，不合并、不修改它们的代码。通用运行时在 OmicsWeave 内独立实现。

## 3. 起点用例：CENOMultiome 与 NonB

### 3.1 两个项目做什么

- **CENOMultiome**：冻结的 DNA encoder 加多组学 any-to-any 模型，DNA + 已观测组学 → 缺失组学。54 个细胞，hg38，131K（以及 512K / 1M）窗口。
- **NonB**：Non-B DNA（G4、R-loop）。分两条管线：peaks 是 5 kb 窗口二分类；tracks 是 131 kb 窗口的 bigWig 信号回归，最多 146 路 track。数据与训练编排由项目自己的数据层完成。

### 3.2 a 端：涉及的输入格式

| 格式 | 内容 | CENO | NonB | 现在的读法 |
|---|---|:---:|:---:|---|
| FASTA + `.fai` | hg38 参考序列 | ✓ | ✓ | CENO 经 fai 整条读出后离线写入 zarr；NonB 按窗口在线读取 |
| bigWig | 1 bp 信号：ATAC、DNase、RNA（分链）、CAGE（分链）、ChIP fold change、G4 / R-loop CUT&Tag | ✓ | ✓ | CENO 离线转成整染色体 zarr（bf16，或 128 bp 均值）；NonB 用 pyBigWig 按窗口在线读，文件在远端对象存储上 |
| BED | peaks 区间集合 | — | ✓ | 内存中合并 / 求交重复、扣背景、二分查询、栅格化成 0/1 bins |
| `.mcool` / `.cool` | Hi-C 接触矩阵（HDF5） | ✓ | — | h5py 读 1 kb 分辨率的 cis 区域，使用 balance 权重 |
| `.hic` | Hi-C 接触矩阵（Juicer） | ✓ | — | 离线构建 targets |
| `.pairs` | Hi-C 接触对 | 元数据层面 | — | 下载时作为候选格式 |
| BAM | 比对读段 | ✓ | — | 离线转 bigWig，**通过子进程调用脚本** |
| CTSS `.bed.gz` | FANTOM CAGE 起始位点计数 | ✓ | — | 离线转成正负链 bigWig |
| tsv / json 清单 | 文件、细胞、track 的元数据与分组 | ✓ | ✓ | 各自的清单格式（CENO `contract.json` 等；NonB `DataPaths` JSON） |

CENO 的自建中间产物：`genome.zarr`（整染色体数组，zstd 分块）、`reference.zarr`（留一 atlas 均值）、`nonzero/*.json`（全基因组非零均值）、Hi-C targets 的 `.npz`。

### 3.3 b 端：训练样本的形态

两个项目的训练样本都以**基因组窗口**为单位，同一窗口上叠加多层坐标对齐的数据：

| 层 | CENO | NonB |
|---|---|---|
| 序列 | uint8 ASCII tokens `[L]` + valid，L = 131072（DNA 上下文可扩到 512K / 1M，监督中央 131K） | peaks：one-hot `(L, 4)`，L = 5000；tracks：字符串，L = 131000，送 CENO tokenizer |
| 1 bp 轨道 | 每个细胞最多 9 路（ATAC、DNase、CAGE±、RNA polyA±、RNA total±），各带 valid mask | `(L, n_tracks)` 的 bp 级信号，n_tracks 最多 146 |
| 粗分辨率轨道 | 128 bp 的 histone（11 种 mark）和 CTCF，各带 valid | 训练端 `avg_pool1d` 成 `resolution` bins |
| 二维 | Hi-C 接触矩阵（2048 bp bins）+ mask；1 kb / 5 kb 输入 | — |
| 标签 / 掩码 | observed / target 掩码（任务格点）、dropout 区间、每 track 缩放均值、留一 atlas 参考 | 窗口标签（任一 track 覆盖即为 1）、`(n_tracks, nbins)` 0/1 覆盖 |

### 3.4 共同的样本构造逻辑

- **窗口采样**：CENO 用固定锚点区间加 shift 增强；NonB tracks 整染色体平铺、相邻窗口随机重叠、每 epoch 用 `seed + epoch` 重新覆盖；NonB peaks 正负例采样，按 locus 分折。划分按整条染色体或 locus 分组，防止泄漏。
- **逐层变换**：归一化（bigWig 头部 `sumData / nBasesCovered`，或全基因组非零均值）、重复聚合（max / mean / merge / intersect）、背景扣除（`log2((a+1)/(b+1))`）、分 bin 求均值、bf16 量化、shift 时 128 bp bin 的线性插值。
- **反向互补**：序列取互补并翻转，同时交换正负链 track。
- **多细胞共享序列**：CENO 中同一位点的多个细胞共用同一段序列（`sequence_key`）。

### 3.5 对 OmicsWeave 的启示

1. **item 是「坐标系上的区间查询」，不是行号。** 访问函数给出 `(contig, start, end, shift, strand, cell/sample)`；不同于 CellPack 每个 cell 对应独立记录，这里相邻窗口、平铺重叠和 shift 会大量共享底层数据，按基因组分块缓存的收益很高。
2. **主流格式自带两级结构，适合「原始块 → 解码块」的缓存模型。** bigWig 是 R-tree 索引加 zlib 压缩数据块，`.hic` 是块索引加压缩块，BGZF / tabix 同理，FASTA 由 fai 直接给出字节偏移。这些格式可以按原样做显式 I/O 和预取，不必先转成自有格式。
3. **CENO 离线转 zarr、NonB 在线读 pyBigWig，是同一个问题的两种妥协。** 前者产生海量分块文件，还要整套清单来管理；后者每个 worker 进程逐窗口读远端存储上的 bigWig，按 float64 逐碱基估算内存，需要大量 CPU 和内存。OmicsWeave 要做的是直接从 bigWig 按已知访问序列预取和缓存，让在线读取够快，从而不需要离线转换。
4. **变换应当是可组合的层操作。** 归一化、重复聚合、背景扣除、分 bin、量化、RC 换链在两个项目里反复出现，属于 M 端（或 a、b 之间）的通用 process，不应由每个项目各写一份。
5. **BAM → bigWig 目前通过子进程完成**，违反 1.2，是用 htslib 进程内实现的候选。
6. **Hi-C 是二维 item**：一个窗口对应一个 `n × n` 的子矩阵，输出形状和缓存粒度都与一维轨道不同，需要单独设计。

## 4. 代码组织

### 4.1 格式组件

- 每种 a 端格式在 `crates/` 下有独立的 crate，目录用短名，包名带 `omicsweave-` 前缀：例如 `crates/bbi` 即 `omicsweave-bbi`。按容器格式命名，而不是按其中一种文件类型：bigWig 和 bigBed 共用 BBI 容器，因此是 `bbi`。
- 代码量小的格式用纯 Rust 重写（例如 bbi）；需要沿用 C 实现的走 `-sys` 惯例：`-sys` crate 只负责编译和链接，安全 API 放在不带 `-sys` 的 crate 中。
- 第三方 C 源码放在使用它的 crate 内部的 `vendor/` 下，按 release 原样拷贝，旁边的 `UPSTREAM.md` 记录上游地址、版本、校验和与许可证。放在 crate 内部，crate 才能独立打包和被引用。压缩库直接依赖 `libz-sys`、`libdeflate-sys` 等现成的 static `-sys` crate，不自行 vendor。
- 组件接口遵循 1.5：解析元数据得到需要读取的范围和需要分配的大小，在调用方的 buffer 上解码。`omicsweave-bbi` 是第一个按此实现的组件。

### 4.2 对齐测试

- 参照实现放在 `refs/<名称>/`，是 `publish = false` 的 `-sys` crate（包名 `<名称>-sys`），只出现在格式组件的 `[dev-dependencies]` 中，不进入 wheel。它只提供三样东西：`vendor/` 中的上游源码、`build.rs` 构建，以及与上游头文件一一对应的原始 FFI 声明。不包含任何测试逻辑或封装。
- 测试逻辑全部写在被测组件的 `tests/` 中：差分测试、测试辅助代码，以及对参照实现原始 FFI 的封装（例如 `crates/bbi/tests/common/libbigwig.rs`）。
- 测试直接在 Rust 里调用双方的函数接口做差分对比（proptest 随机生成输入和查询），不通过 CLI，也不比较输出文件。
- 参照实现选择现有 Python 管线实际使用的库，以保证替换后行为不变：bbi 对齐 libBigWig（pyBigWig 的底层）。
- 参照实现无法生成的输入（例如 libBigWig 不能写 bigBed），由测试代码中的最小写入器生成，再交给参照实现和我们的实现分别读取。参照实现是独立的读取器，两边一致即同时验证了写入器和读取器。
- 参照实现自身的已知缺陷不作为对齐目标。例如 libBigWig 写入 zoom 记录时，缓冲区写满会在保存前重置累计和，导致记录的 `sum` 不精确；zoom 测试因此只比对覆盖碱基数、min、max。

### 4.3 标准访问与输出（`omicsweave-core`）

同一系列的不同格式，访问方式和输出布局完全相同；格式之间的差异只存在于各自的低层步骤（算字节范围、解码块）中，这些步骤保留给自行调度 I/O 的 loader 使用。

**访问 key**：`Region { chrom: &str, start: u32, end: u32 }`，0-based、左闭右开。坐标统一为 `u32`，与 BBI、2bit、bigBed 等格式的存储宽度一致；文件偏移等元数据仍为 `u64`。

**四个系列**：

| 系列 | 格式 | 访问 trait | 输出 |
|---|---|---|---|
| Sequence | `.fa`、`.fa.gz`、`.2bit` | `SequenceRead` | `[u8]` ASCII，长度 = `region.len()` |
| Track | bigWig | `TrackRead` | values：`[f32]`，长度 = `region.len()`，无数据为 NaN；entries：`[TrackEntry]`，**裁剪到查询范围** |
| Interval | bigBed、BED | `IntervalRead` | `[IntervalEntry]` + 注释字符串列；**不裁剪** |
| Contact | `.hic` | `ContactRead` | 稠密 tile：行优先 `[f32]`，形状 = 行 bin 数 × 列 bin 数；或 `[ContactEntry]` 非零格（按 (row, col) 排序） |

- **行优先的定长 entry**：`TrackEntry { start: u32, end: u32, value: f32 }`、`IntervalEntry { start: u32, end: u32, annotation: u32 }`，都是 12 字节的 `#[repr(C)]` 结构，与 numpy 结构化 dtype 一一对应。
- **注释字符串列**：Arrow `LargeString` 布局（`i64` offsets 长度 n+1，加 `u8` data）。`IntervalEntry::annotation` 是该列中的行号，所以对 entry 排序时不必移动字符串。每行内容是记录第三列之后的全部字段，用 tab 连接，与 bigBed 的存储方式相同。
- **Track 裁剪、Interval 不裁剪**：Track 的区间只是信号的一种表达方式，不应影响数据内容，因此查多少给多少；Interval 的每条记录是完整的特征（例如一个 peak），注释描述的是整个区间，因此保持原样。
- **变长输出：一次解码，暂存，再拷贝**：`TrackRead::query_entries` / `IntervalRead::query` 对查询区间只解码一次，把结果暂存在 reader 内部并返回输出大小；调用方分配好内存后，`write_entries` / `write`（`&self`，不再读取或解码）把暂存结果拷贝进调用方的缓冲区。暂存属于 reader 的内部工作状态，交付的输出始终在调用方的内存中。缓冲区大小不符返回 `Error::OutputSize`，未查询就写返回 `Error::NotQueried`。
  - bigBed：查询时把解压后的块依次追加到一块只增不缩的内部缓冲区，暂存的记录只保存注释在其中的位置，不拷贝注释字节；`write` 是唯一一次拷贝。
  - bigWig entries：暂存解码并裁剪好的 `TrackEntry`，`write_entries` 是一次 memcpy。
  - BED：注释常驻在表中，查询只暂存命中记录的下标，`write` 不再重做查找和筛选。
- **越界与未知染色体**：Sequence 的输出必须定长，区间超出序列即报 `RegionOutOfRange`；Track 和 Interval 允许超出染色体末端，没有数据的部分就是没有。未知染色体：Sequence 和 Track 报 `UnknownChrom`；Interval 返回空，因为 BED 没有染色体列表，bigBed 与之保持一致。空区间在所有系列中都返回空。
- **统一 Error**：所有格式 crate 使用 `omicsweave_core::Error`，以后 Python 绑定只需映射一套异常。
- **字节源**：读取器通过调用方提供的 `ReadAt`（`[u8]`、`Vec<u8>`、`FileSource` 等）读取所需字节；I/O 来自何处、是否缓存，由调用方决定。

**不重复解码**：任何重复扫描或重复解码，只有在开销小于一次拷贝时才允许；否则一律解码一次、暂存、再拷贝。定长输出（Sequence、Track values）直接解码进调用方的缓冲区，不经过暂存。

**Contact（Hi-C）**：只支持 `.hic` v8 和 v9，其它版本报 `Error::Unsupported`。只提供原始接触计数：归一化和 expected 属于分析逻辑，格式层不读取，也不计算。查询是 `ContactQuery { rows, cols, resolution, unit }`。

- **区间与方向**：行和列可以在不同染色体上，区间覆盖的每个 bin 都计入 tile。无论文件存的是对称矩阵的哪一半，都按查询方向给出完整的 tile（对角线另一侧自动镜像）。
- **单位 `unit`**：`Bp` 时区间是碱基坐标；`Frag` 时区间是酶切片段编号，片段 `i` 从第 `i-1` 个酶切位点到第 `i` 个。
- **稠密输出**：没有接触的格为 0。
- **entries**：只含计数有限且非零的格，坐标为 bin 起点（单位同查询）。
- **只解码一次**：块索引和酶切位点首次使用时读入并缓存；一次查询中每个块只解压、解码一次。
- **v8 与 v9 的差别**只在四处：header 中的 NVI 指针和染色体长度宽度、footer 中 `nBytesV5` 的宽度、块内坐标宽度标志（v9 才有）、染色体内的块编号（v9 沿对角线按深度编号，受 `v9-depth-base` 控制）。
- 读出的计数以 straw 的 `observed` / `NONE` 为参照逐位校验。

**写 `.hic`**：
- `HicWriter` 按 Juicer 的布局写 v8 或 v9（默认 v9）：header、blocks、矩阵记录、footer、normalized expected、norm vector index。expected 和归一化列表写为空，straw 仍可用 `NONE` 读取。
  - 调用方给出分好 bin 的矩阵。
  - 每个 block 选用最小的表示：list-of-rows 或 dense，计数用 short 或 float；v9 中坐标可在 short / int 之间按需选择。
  - 支持 v9 的 `v9-depth-base` 布局。
  - **块大小可配置**（`BlockLayout` / `BlockBins`，`HicWriter::write_matrix_with_blocks` 可按单个矩阵覆盖）。默认染色体内矩阵用 128 bin 的块，染色体间矩阵沿用 Juicer 网格（约 1000 bin）。理由：按接触密度 ∝ 1/距离 估算，1 kb 下一个 131 kb 的 cis 窗口在 v9 对角线布局中，用 1000 bin 块要解码约 15 倍于所用的接触，用 128 bin 块约 2.5 倍、只涉及约 2 个块；更小的块主要增加请求数。染色体间矩阵稀疏，块切小只会产生大量碎块。
  - 块大小写在每个矩阵记录里，任何读取器都按文件中的值解析；straw 对 1–400 bin 的块与我们逐位一致。注意 straw 计算 v9 染色体内块号时固定用以 2 为底的 log 深度，忽略 `v9-depth-base`，所以非默认深度布局的文件只有我们的读取器能正确读。
- `HicBuilder` 从接触对构建完整文件，相当于 Juicer 的 `pre`：生成所有 BP / FRAG 分辨率，不计算 expected 和归一化。
- 写出的文件以 straw 为参照逐位校验。

**个体化序列（consensus）**：`omicsweave-consensus` 是 `bcftools consensus` 1.23 的移植，源自 pyconsensus。
- **输入与输出**：VCF/BCF 由 `omicsweave-vcf` 一次读入，常驻内存（`VcfStore`，只读、可跨线程共享）。每条记录只留构建要用的部分：位置、参考跨度、等位，以及按 bcftools 读法保存的 `int8` 基因型；常见的双等位基因型压成每个值 2 bit。参考序列经任意 `SequenceRead` 读取。输出长度随 indel 变化，沿用"查询返回大小、再拷进调用方缓冲区"的模式。
- **批量**：`load` 读一次区间的参考序列（含 mask），并找出该区间的记录；之后每次 `query` 按 `Alleles` 生成一条序列（对应 `-s`、`-H`、`-I`）。同一区间的大量样本、单倍型共用这一次加载，每条序列的代价只剩线性拷贝和逐记录的等位选择。
- **单一执行路径**：逐条移植 `apply_variant` 的状态机，不设快速路径分支。consensus.c 原地修改缓冲区，每遇到 indel 就整体平移尾部；这里改为"已改写的前缀 + 未改动的参考后缀"，记录按序到达，且不会改写前一条记录之前的位置，所以每个区间的耗时与长度成线性。
- **支持的选项**：`bcftools consensus` 的选项全部支持，包括 `-I`、`-M`、`-a`、`-m`/`--mask-with`、`--mark-*`、`-c`、`--regions-overlap` 和 `-i`/`-e`。
- **过滤**：consensus.c 对不满足 `-i`/`-e` 的记录只是跳过，所以过滤在 `VcfStore` 读入时一次做完（`omicsweave-vcf` 的 `Filter`，即 `filter.c` 的重写）。表达式编译失败时读入报错；某条记录求值出错，则记在该记录上，查询经过它时才报，与 bcftools 读到该记录才报错的时机一致。
- **不发明缓存格式**：不使用 pyconsensus 的 `.cvcf`。只解析一次就能服务所有构建；多进程场景在 fork 前加载，需要更快的加载时用 BCF。
- **对齐基准**：`refs/bcftools` 在进程内运行原版 bcftools，以它的输出为准逐字节对齐，包括它的若干行为怪癖。

**BED 内存对象**：BED 不支持随机访问，一次性载入（含 gzip 流式解压）时直接构建拆分格式：entry 按染色体分组、组内按 start 稳定排序，附带前缀最大 end 用于二分查找重叠区间，注释按文件顺序存为字符串列。查询结果与 bigBed 完全相同。

### 4.4 Sidecar（`<源文件>.sidecar`）

loader 要在读取数据块之前静态分配内存，因此需要知道每个中间产物的大小。有些格式没有记录块解压后的大小：bigWig / bigBed 的 R-tree 和 `.hic` 的块索引只给出压缩字节数，zlib 流本身也不记录原始长度；`.hic` 的块编码（list-of-rows 或 dense，short 或 int）还是按块选择的。sidecar 为每个块记录这些信息，属于 1.1 所说的派生产物：可以随时删除，并按源文件身份重建。

- **格式**：紧凑的定长小端二进制，不压缩。`omicsweave-core::sidecar` 定义公共的 48 字节容器头：magic `OWSIDECR`、容器版本、格式种类、body 版本、源文件身份、body 长度和校验和。body 由各格式定义，由 8 字节对齐的定长记录表组成，读入后可以直接按下标访问。
- **源文件身份**：文件长度，加上对开头 64 KiB、结尾 64 KiB 和 14 个均匀分布的 4 KiB 页做的指纹。它只依赖文件内容，不依赖路径或 mtime，所以源文件和 sidecar 一起复制、移动后仍然有效。身份不符时报 `Error::StaleSidecar`，调用方应据此重建。
- **bbi**（`BbiSidecar`）：base 级别和每个 zoom 级别的所有块，每块记录 span、文件偏移、压缩大小、解压大小和条目数。块在级别内的位置就是它的稠密序号。`BigWigWriter` / `BigBedWriter` 在写文件时同时产出 sidecar，`BbiSidecar::build` 从已有文件构建。bigWig 也可以不用 sidecar，按头部的 `uncompressBufSize` 上界分配。
- **hic**（`HicSidecar`）：每个 (chr1, chr2, unit, bin size) 矩阵的所有块，每块记录块号、位置、压缩大小、解压大小和解码出的记录数。`HicWriter` / `HicBuilder` 的 `finish_with_sidecar` 在写文件时同时产出 sidecar，`HicSidecar::build` / `build_filtered` 从已有文件构建（后者只处理选中的矩阵）。

### 4.5 FASTA 索引（`.fai`、`.gzi`）

`.fa` / `.fa.gz` 的索引沿用标准格式，与 `samtools faidx` / `bgzip` 逐字节一致，接口与 sidecar 同一风格：`Fai` / `Gzi` 各有 `build`（从源文件构建）、`parse` / `to_bytes`、`save` / `load`（`<源文件>.fai`、`<源文件>.gzi`，原子写入）和 `build_file`。

- `build_index` 对 bgzip 文件只解压一遍，同时得到 `.fai` 和 `.gzi`；`Gzi::build` 只读每个块的头部和 ISIZE，不解压。
- `FastaWriter` 写文件时直接按写出的布局记下 `.fai` 条目（压缩时同时记下 `.gzi`），不重新解析写出的内容。
- 空的 BGZF 块（包括 EOF 块）不进 `.gzi`，与 htslib 一致。

### 4.6 格式转换（`omicsweave-convert`）

`omicsweave-convert` 是端到端的转换工具：在给定的**线程数和内存量**内，流式、并行地把 a 格式文件写成 b 格式文件。

| 转换 | 输入 | 输出 | 对齐基准 |
|---|---|---|---|
| `FastaToBgzipFasta` | `.fa`（纯文本 / gzip / BGZF） | `.fa.gz` + `.fai` + `.gzi` | 与 `compress_fasta`（即 bgzip 的块划分）逐字节一致 |
| `FastaToTwoBit` | `.fa` | `.2bit` | 与 `faToTwoBit` 布局逐字节一致，lib2bit 读回 |
| `VcfToBcf` | `.vcf` | `.bcf`（+ `.csi`） | 未压缩流与 htslib `wbu` 逐字节一致；索引与 `bcftools index` 同一内容 |
| `SamToBam` | `.sam` | `.bam`（+ `.bai` / `.csi`） | 未压缩流与 htslib 逐字节一致；索引与 `samtools index` 同一内容 |
| `BedToBigBed` | `.bed` | `.bigbed` + `.sidecar` | 与 `BigBedWriter` 逐字节一致 |
| `CoolToHic` | `.cool` / `.mcool` | `.hic` + `.sidecar` | 与 `HicWriter::write_matrix` 逐字节一致 |

**选择转换。** `Format` 列出 OmicsWeave 读写的全部格式（含 `.fai`、`.gzi`、`.csi`、`.bai`、`.sidecar`），同一格式的不同存储（plain / gzip / bgzip）各是一个值，名称与 Python 的 `File.format` 一致；`Format::detect` 按内容、`Format::from_extension` 按名称识别。`Conversion::between(from, to, src, dst)` 按格式对选择转换，源格式的任一存储都可以；`Conversion::infer` 两端都按扩展名判断；`Conversion::configure(&Options)` 应用选项，不适用于该转换的选项报错。`Report.files` 带每个文件的 `Format`。

**分工。** Python（`omicsweave.convert`、`ow convert`）只做顶层组装：把路径或文件对象交给原生层，把返回的格式名组装成文件对象，命令行负责参数、进度与输出。识别格式（源按内容、目标按名称）、校验选项、预算缺省值、解析大小都在 Rust。

**代码结构。** `format`、`conversion`（选择、`plan` / `run`、运行上下文）之外，每种源格式一个模块（`fasta`、`vcf`、`sam`、`bed`、`cool`），共用部件在 `pipeline/`：`input`（压缩识别、并行解压、按行分批）、`output`（临时文件与原子改名、spill 目录）、`sort`（外排序）、`records`（VCF → BCF、SAM → BAM 共用的「文本行 → BGZF 记录 + 分箱索引」流水线）。

**三段式流水线。** 每条转换都由同样的三种阶段组成：

1. **串行读取 / 分帧**：把输入切成单元（行批、BGZF 块、像素区间），并记下只依赖未压缩内容的元数据（`.fai`、zoom 汇总、块号）。开销小。
2. **并行纯计算**：解压、解析、编码、压缩，以任务的形式提交到共享的 `omicsweave_core::exec::Pool`。任务只拥有自己捕获的数据（`'static`），不借用 writer 的状态。
3. **串行按序提交**：按提交顺序取回结果，写出，得到压缩后的偏移，建立依赖偏移的索引（`.gzi`、R-tree、hic 块索引、BAI/CSI、sidecar）。

`core::exec` 的 `Ordered` 把这一接缝标准化：`Inline` 在提交时立即执行，`PoolQueue` 在线程池上执行，结果都按提交顺序返回。各格式的 writer 只保留一条代码路径，`Inline` 只是执行器的一种。由此得到**确定性不变式**：输出与线程数无关，也与排序是否 spill 到磁盘无关，测试对每条转换都检查这一点。

**哪些能并行，哪些不能。**

| 转换 | 并行 | 串行 |
|---|---|---|
| FASTA → bgzip | deflate；BGZF 输入的解压 | 读取、`.fai` 扫描、写出 |
| FASTA → 2bit | 各片段打包，统计 N / mask 游程 | 解析 FASTA、合并游程、目录 |
| VCF → BCF、SAM → BAM | 按行批解析编码、deflate；BGZF 输入的解压 | 切行、BGZF 分帧、建索引、写出 |
| BED → bigBed | 解析、分区排序、块压缩 | 染色体编号、zoom sweep、R-tree |
| cool → hic | chunk 读取与解码、像素定位、分区排序、块编码与压缩 | 矩阵记录、footer |

普通 gzip（非 BGZF）无法拆分，只能在单独的线程里流水线解压，不能并行。

**BGZF（`omicsweave-bgzf`）。** fasta、vcf、sam 共用。`Writer` 按 htslib `bgzf_write` / `bgzf_flush_try` / `bgzf_flush` 的规则切块，所以未压缩块与 htslib 一致。块在执行器上压缩；写入时报告的位置是 `Logical` 偏移（块序号代替压缩地址），文件写完后由 `Blocks::resolve` 换算成虚拟偏移。`BinningIndex` 按 `hts_idx_push` / `hts_idx_finish` 的逻辑建 BAI / CSI：先在 `Logical` 偏移上做在线合并，偏移换算之后再执行那些需要比较地址大小的步骤（合并小 bin、合并同块 chunk）。`Reader` 在执行器上预读并解压后续块。

**预算。** `Budget { threads, memory, tmp_dir }`。内存采用规划式分配，运行时不计量：`plan` 在读数据之前给出 `MemoryPlan`，包括三部分：
- `fixed`：常驻部分，如索引、header、bins 表；
- `pipeline`：每个队列的在途上限 × 单元大小；
- `buffer`：需要重排的转换使用的排序缓冲区。

`fixed + pipeline` 超出预算时直接报错，不在运行中途 OOM。批大小随预算自适应（64 KiB 到 4 MiB）。核算只覆盖大块缓冲，分配器开销、线程栈和小元数据不计入，RSS 允许超出约 10% 加 64 MiB。

**超过内存时。**
- **纯流式**（FASTA → bgzip、VCF → BCF、SAM → BAM）：工作集只有在途批次，任意大小都能处理。
- **单元级**（FASTA → 2bit）：2bit 记录头在打包数据之前，所以一条序列的打包数据先暂存；超出缓冲时写到临时文件。目录在所有记录之前，序列名来自源文件的 `.fai`，没有 `.fai` 时多扫一遍输入。
- **需要重排**（BED → bigBed、cool → hic）：用通用外排序 `pipeline::sort::Sorter`。缓冲满时，并行排序各部分，归并写成 run 文件，最后与剩余的内存部分做 k 路归并。键唯一（BED 加输入序号），所以结果与在哪里 spill 无关。cool → hic 按 (染色体对, 块号, y, x) 排序，之后用 hic 的流式接口（`begin_matrix` / `add_block` / `end_matrix`）逐块写出，整张矩阵不需要放进内存。

spill 文件默认放在输出文件旁边（不用 `/tmp`）。所有输出先写到 `.partial`，整个转换成功后才改名；失败时清理。

**内存数据结构。** 每种格式只定义转换需要的流式单元，不定义整文件模型，也不对外提供操作接口：

| 数据模型 | 单元 |
|---|---|
| VCF / BCF | 复用 `Record` |
| SAM / BAM | 行到 BAM 记录的编码器（`BamEncoder`） |
| FASTA / 2bit | `TwoBitPiece` / `TwoBitRecord` |
| Interval | `(chrom, start, end, rest)` |
| Contact | 按 `MatrixLayout` 定位后的 `Record` |

查询接口仍然只有"a 文件 → b 内存产物"一种；Python 层只暴露 `omicsweave.convert`（`run` / `plan` / `routes`），命令行是 `ow convert`。

**与参照工具的有意偏差。**
- VCF 中 header 未声明的名字报错。htslib 会动态补 header，但写出的 BCF 会损坏；需要补 contig 时开启 `add_missing_contigs`，先扫一遍 CHROM 列。
- SAM 中 RNAME 不在 @SQ 里时报错。htslib 只警告，并把读段当作未比对。
- BED 的染色体按首次出现顺序编号；bigBed 的 B+ 树按名字索引，读取不受影响。
- cool 只取上三角（`square` 模式同样如此），不读 balance weight。
- `.hic` 不加 Juicer 的 `All` 染色体。

### 4.7 Python 层（`omicsweave`）

Python 绑定端到端处理加载：按路径打开文件、按 key 查询、拿到 numpy 结果，不暴露索引、数据块、解压器等中间数据。

- **类型都定义在 Python 里**：key（`Region`、`ContactWindow`、consensus 的 `Haplotype` / `Samples` / `AllAlt`）、结果（每个结果带上它的 key，数据全是 numpy 数组）、选项和异常都是 Python dataclass 或异常类；原生模块 `omicsweave._native` 只收发基本类型和 numpy 数组。
- **三层**：`types` 是与格式无关的 key 与结果；`formats` 每种文件格式一个类，文件对象只描述路径上的格式，同一格式的不同存储（plain / gzip / bgzip）各是一个类，索引与伴随文件（`.fai`、`.gzi`、`.sidecar`）也是平等的格式类，writer 属于各自的格式；`query` 是读取器，接收文件对象，按 key 返回 `types` 的结果。格式按内容判断（原生 `detect`，gzip 文件解压开头几个字节），与扩展名无关。
- **列优先**：区间、track 区间、接触 cell 按列输出成多个连续数组，不转成 Python 对象；字符串列用 Arrow `LargeString` 布局（`offsets: int64`、`data: uint8`），可零拷贝转成 pyarrow。core 的行优先 entry 只在原生层内部暂存，拆列是那一次拷贝。
- **内存**：定长结果可由调用方传入 `out=`，原样写入；变长结果在原生层按查询得到的大小分配，所有权直接交给 numpy，不再拷贝。
- **伴随文件只显式构建**：FASTA 的 `.fai`（bgzip 还有 `.gzi`）是按 key 查询的前提，可显式给出路径，否则在源文件同目录按同名查找；打开时从不自动构建，缺失即报错，由 FASTA 文件对象的 `build_index()` 手动构建或由 `FastaWriter` 写出。按 key 直接查询不读 sidecar：从已有文件建 sidecar 要遍历整个文件，单次查询用不上；它只供以后的 loader 使用，由文件对象的 `build_sidecar()` 或 writer 的 `sidecar=True` 显式产出。
- **并发**：每个读取器内部状态由一把锁保护，只在释放 GIL 后加锁，解码期间不持有 GIL；读取器可 pickle，按文件和选项重新打开。
- **异常**：`Error` 的每个变体对应 `omicsweave.errors` 中的一个异常类，字段变成属性，并同时继承语义最近的内置异常。

## 5. 待设计

- Source 接口：catalog、对象依赖、decode kernel 的 Rust trait 与 C ABI
- Planner 的通用对象模型与代价模型
- Output schema：规范行的类型体系与变长输出
- Python 层之上：数据源声明（清单 / 配置）、访问函数构造、与 PyTorch 的衔接
- 派生产物的缓存管理与失效策略

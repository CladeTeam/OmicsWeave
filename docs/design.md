# OmicsWeave 设计

OmicsWeave 是一个系统性的 a2b 框架：a 端是组学数据的输入格式，b 端是模型训练使用的数据形式。框架负责把任意受支持的 a 端数据组织成 b 端训练数据，并把多种组学数据整合到统一的加载器中。

## 1. 设计原则

### 1.1 沿用现有格式，只在真实瓶颈处设计新格式

a 端格式按原样读取，不要求用户先把数据转换成 OmicsWeave 自有的存储格式。很多生信格式效率不高，例如普通 gzip 压缩的 FASTA 无法随机访问，压缩收益也有限，但只是格式不理想不构成重新设计的理由。

只有当某个 a2b 路径在真实负载下出现可测量的瓶颈，且索引、缓存等手段无法解决时，才设计新的存储格式。CellPack 就是这样产生的：单细胞表达矩阵的 shuffle 训练加载被 H5AD 的读取方式卡住，才有了逐细胞可独立解码的新格式。

由此得到以下约束：

- **原始文件始终是事实来源。** 框架不改写、不替换用户的原始数据。
- **允许派生产物，但它们只是加速手段。** 索引（`.fai`、`.tbi`、`.csi`）、解压副本、编译缓存（例如 pyconsensus 的 `.cvcf`）都属于派生产物：可以随时删除、按需重建，并按源文件身份（路径、大小、mtime 或内容哈希）校验是否失效。派生产物不得成为唯一的数据来源。
- **查询只覆盖能高效随机访问的格式，其余格式只做转换源。** 按区间查询的读取器只读带索引、可按块定位的格式：bgzip FASTA / 2bit、bigWig、bigBed、`.hic`、BAM（`.bai` / `.csi`）、CRAM（`.crai`）、BCF（`.csi`）。BED / BED.gz 这类直接存文本、只能整体扫描的格式，以及只能整块读取、不能只读一个 tile 的 cooler，不提供查询路径，由 `ow convert` 转成上述格式（BED → bigBed、`.cool` / `.mcool` → `.hic`、VCF → BCF、SAM → BAM）再读。
- **转换只做格式之间的变换。** 覆盖度计数、归一化、过滤等属于业务逻辑，不进入转换：例如 BAM 到 bigWig 的覆盖度不是格式转换，由调用方自己计算。
- **低效格式的处理顺序**：先转换成带标准索引的格式；仍有可测量的瓶颈时，再构建派生缓存；仍不够时，才进入新格式设计。
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

### 1.5 Rust 底层组件在调用方分配的 buffer 上操作

Rust 底层格式组件的输入和输出使用调用方分配的 buffer，组件在上面读写，由上层运行时或绑定层管理交付数据的生命周期。组件保留外部分配接口，供 loader 的 cache frame、输出 ring 等受预算管理的内存使用。

Python 查询层是这些组件的调用方：它负责分配输出、管理 I/O 和解码工作区，并一次返回结果。Python 用户不需要先查询输出大小、申请 buffer，再执行写入。定长查询的 `out=` 只是已有批次内存需要复用时的可选入口，见 4.7。

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

**访问 key**：坐标统一为 `u32`，0-based、左闭右开，与 BBI、2bit、bigBed 等格式的存储宽度一致；文件偏移等元数据仍为 `u64`。key 分两种：

- `Region { chrom, start, end }`：一条染色体上的一个区间，是稀疏查询（特征、读段、记录、接触像素）的 key。
- `Locus { chrom, spans, strand }`：同一条染色体上升序、互不重叠（可以相接）的若干片段，加一条链 `Strand`（`+`、`-`、`.`、`?`），是稠密查询（碱基、逐位信号、bin 统计、接触 tile）的 key。`Region` 可直接转成单片段、无链的 locus；剪接后的转录本就是由外显子组成的 locus。

**方向与拼接的统一规则**：

- **稠密输出**按 locus 的 5'→3' 写出：各片段首尾拼接，`-` 链整体翻转，碱基再取互补（IUPAC 码映射到互补集合的码，保留大小写，`U` 按 `T` 处理）；`.` 和 `?` 按 `+` 读。长度 = 各片段长度之和。
- **稀疏输出**带自己的坐标：保持基因组坐标和文件顺序，strand 既不过滤，也不翻转。
- **只解码一次**：先求出与任一片段重叠的数据块并去重，每块只读、只解压一次，再分发到拼接后的输出位置。不会去解码片段之间的内含子；多个片段共享同一个块时也不重复解码。bgzip FASTA 中，共享 BGZF 块的相邻片段合并成一段一起解码。
- **翻转在原地进行**：先按正向解码进调用方的 buffer，再原地 reverse（`core::orient`），不需要额外内存。

**五个系列**（另有 BCF 的基因型矩阵，见下文「变异（BCF）」）：

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
- **统一 Error**：所有格式 crate 使用 `omicsweave_core::Error`，以后 Python 绑定只需映射一套异常。
- **字节源**：读取器通过调用方提供的 `ReadAt`（`[u8]`、`Vec<u8>`、`FileSource` 等）读取所需字节；I/O 来自何处、是否缓存，由调用方决定。

**不重复解码**：任何重复扫描或重复解码，只有在开销小于一次拷贝时才允许；否则一律解码一次、暂存、再拷贝。定长输出（Sequence、Track values）直接解码进调用方的缓冲区，不经过暂存。

**bin 汇总（bigWig）**：`TrackRead::bins(locus, BinSpec { resolution, phase, stats, source }, out)` 在读取器内把信号汇总成定宽 bin，不展开逐碱基数组。

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
- **输入与输出**：只读 BGZF BCF 与它的 CSI，经 `omicsweave-vcf` 的 `BcfReader` 按区间读取（见下文「变异（BCF）」）；VCF 文本先 `ow convert x.vcf.gz x.bcf --index`。`load` 只解码索引给出、且与区间重叠的记录，每条一次；每条记录只留构建要用的部分：位置、参考跨度、等位，以及按 bcftools 读法保存的 `int8` 基因型。取哪些记录与 htslib synced reader 对 `--regions-overlap` 的判断一致：先按 `bcf_itr_queryi` 的存储跨度取候选，再按解码出的跨度（`0` 看起点，`1` 看整条记录，`2` 看去掉与所有 ALT 共同前缀后的部分）筛选，所以起点在区间之前的删除等记录也会被取到。参考序列经任意 `SequenceRead` 读取（`Locus::from(region)`）。输出长度随 indel 变化，沿用"查询返回大小、再拷进调用方缓冲区"的模式。
- **链**：`write(strand, out)` 按链给出序列：`-` 时边拷贝边反向互补（IUPAC 码取互补集合，保留大小写），其余链按正向输出；chain 始终是正链坐标。Python 层按 `Region.strand` 传入，批量结果中每条序列在自己的槽位内反向；`Spliced` key 不适用（consensus 只编辑一段连续区间），报 `TypeError`。
- **批量**：`load` 读一次区间的参考序列（含 mask），并解码该区间的记录；之后每次 `query` 按 `Alleles` 生成一条序列（对应 `-s`、`-H`、`-I`）。同一区间的大量样本、单倍型共用这一次加载，每条序列的代价只剩线性拷贝和逐记录的等位选择。
- **单一执行路径**：逐条移植 `apply_variant` 的状态机，不设快速路径分支。consensus.c 原地修改缓冲区，每遇到 indel 就整体平移尾部；这里改为"已改写的前缀 + 未改动的参考后缀"，记录按序到达，且不会改写前一条记录之前的位置，所以每个区间的耗时与长度成线性。
- **支持的选项**：`bcftools consensus` 的选项全部支持，包括 `-I`、`-M`、`-a`、`-m`/`--mask-with`、`--mark-*`、`-c`、`--regions-overlap` 和 `-i`/`-e`。
- **过滤**：consensus.c 对不满足 `-i`/`-e` 的记录只是跳过，所以过滤在 `load` 解码记录时做（`omicsweave-vcf` 的 `Filter`，即 `filter.c` 的重写），按文件顺序逐条求值；表达式在构建 `Consensus` 时编译一次，每个区间从未求值过的状态开始，与只含该区间的一次 `bcftools consensus` 一致，结果也不依赖查询顺序。表达式编译失败时构建报错；某条记录求值出错，则记在该记录上，查询经过它时才报，与 bcftools 读到该记录才报错的时机一致。
- **不发明缓存格式**：不使用 pyconsensus 的 `.cvcf`，也不把整个文件读成常驻对象：BCF + CSI 已能按区间只读所需的块。克隆（`clone_context` / `try_clone`）共享 header、索引、编译好的过滤式和选项，各自持有读取与解码状态。
- **对齐基准**：`refs/bcftools` 在进程内运行原版 bcftools，以它的输出为准逐字节对齐，包括它的若干行为怪癖。

**比对（BAM）**：`omicsweave-bam`，以 htslib 的读法为准。
- **打开**：读 header，记录的参考编号以二进制参考列表为准；参考名必须是 UTF-8，重名时按名字找到第一个。索引（BAI 或 CSI）由调用方给出；Python 层按 htslib 的顺序在旁边查找（`<bam>.csi`、`<name>.csi`、`<bam>.bai`、`<name>.bai`），找不到即报 `IoError`，打开时从不构建。
- **查询**：按 `hts_itr_query` 由索引算出 chunk（linear index 与 `loff` 给出 `min_off`，右侧 bin 给出 `max_off`，再去包含、去重叠、合并同块），按 `hts_itr_next` 的规则读：每个 chunk 一次 pread，记录的 `[pos, bam_endpos)` 与区间重叠才返回，遇到其它参考或 `pos ≥ end` 即停。结果与 `sam_itr_queryi` 逐条、逐字段一致，包括放在参考上的 unmapped 记录。随机读不依赖 BGZF 之外的派生产物。
- **记录**：校验与 `bam_read1` / `sam_read1` 相同（长度、参考编号范围、CIGAR 与序列长度）；超过 65535 个操作的 CIGAR 放在 `CG:B:I` 标签里时按 `bam_tag2cigar` 还原，aux 输出不含该标签。
- **输出**：`AlignmentEntry`（32 字节：`start` / `end` 为 `i64`，`end` 即 `bam_endpos`；`mate_start`、`mate_tid`、`tlen` 为 `i32`；`flag`、`mapq`），加上按 `AlignmentFields` 选取的变长列（name、ASCII 序列、quality、`u32` CIGAR、原始编码的 aux），都是 offsets + data。查询时只把选中字段的字节暂存一次，`write` 解包碱基并拷出。
- **建索引**：`build_index` 顺序读整个文件（解压并行），按 `bgzf_tell` 的语义取每条记录的结束位置（记录正好在块尾结束时，位置指向下一个物理块，空块也算），交给与 `SamToBam` 共用的 `BinningIndex`，结果与 `samtools index` 内容相同。
- **对齐基准**：htslib 的 `sam_itr_queryi` 与 `sam_index_build3`（`refs/htslib/src/bam.c` 导出 `bam1_t` 字段和 `sam.h` 中的内联迭代器），覆盖 BAI、不同 `min_shift` 的 CSI、多种压缩级别、文件中间的空 BGZF 块和长 CIGAR。
- **不做**：SAM 不支持随机访问，只作为转换源；不做覆盖度、pileup 或从读段生成 consensus：计数、变异检测属于业务逻辑。

**比对（CRAM）**：`omicsweave-cram`，以 htslib 1.25 的读法为准，经 `AlignmentRead` 给出与 BAM 完全相同的输出。
- **版本与解码**：只读 CRAM 3.0 与 3.1，其它版本报 `Error::Unsupported`。容器、slice、压缩头和各数据序列的编码（EXTERNAL、HUFFMAN、BETA、SUBEXP、GAMMA、BYTE_ARRAY_LEN、BYTE_ARRAY_STOP）用 Rust 解码；块的 gzip 用 libdeflate，bzip2、lzma 用静态链接的 `bzip2-sys`、`lzma-sys`，rANS 4x8 / Nx16、自适应算术编码、fqzcomp 与名字分词器来自 `htscodecs-sys`（htscodecs 1.6.7 原样 vendor，SIMD 解码器按 CPU 运行时选择）。块的 CRC32 在解压时校验，从不解压的块不校验，与 htslib 相同。
- **参考序列**：调用方给出任意 `SequenceRead`（bgzip FASTA、2bit 等），按 `@SQ` 名字取 slice 的参考区间并转大写；与 htslib 一样校验 slice 的参考 MD5（可与 CRC32 一并关闭）。内嵌参考与无参考（`RR=0`）的文件不需要参考；多参考 slice 只读取查询所在参考的序列。
- **索引**：`.crai`（gzip 文本）按 `cram_index_load` 载入成嵌套包含表，起点由 `cram_index_query` 的二分查找给出；Python 层在旁边找 `<cram>.crai`、`<name>.crai`，缺失即报 `IoError`，打开时从不构建。`build_index` 只解码多参考 slice（按记录求各参考的范围，同 `cram_index_build_multiref`），其余条目直接取自 slice 头，解压后与 `samtools index` 逐字节一致。
- **查询**：从索引给出的容器顺序读，按 `cram_first_slice` / `cram_next_slice` 跳过或终止容器与 slice，再按 `cram_get_seq` 用记录的 1-based `apos` / `aend` 过滤；结果与 `sam_itr_queryi` 逐条、逐字段一致，包括多参考 slice、放在参考上的 unmapped、超长 CIGAR 与 `CG` 标签。空区间返回空：htslib 对 CRAM 会返回覆盖该位置的记录，这里与 BAM 和 `AlignmentRead` 的约定一致。
- **记录**：按 `cram_decode_slice` / `cram_decode_seq` / `cram_to_bam` 还原：slice 内 mate 的交叉引用与 TLEN 推算，没有存名字的读段按文件名和记录号命名，读组追加为 `RG:Z`，未存储的 MD / NM 按参考生成（htslib 默认如此），`CG:B:I` 按 `bam_tag2cigar` 移回 CIGAR。碱基经 4-bit 编码往返，与 BAM 读出的一致。
- **按字段解码**：CRAM 的各数据序列分块存放。查询由所需字段推出数据序列集合，再闭包到与之共享块（或 CORE 比特流）的序列与标签，只解压用到的块。只要定长字段与 CIGAR 时，解码 BF、CF、RI、RL、AP、MQ、mate 相关序列与读段特征的位置和长度，由特征推出 CIGAR 和比对终点，不重建碱基，也不需要参考；碱基、质量值、名字与标签所在的块都不解压。htslib 完整解码时，持有参考的 slice 中越过参考末端的读段终点算法不同；这里按"完整解码是否持有参考"保留这一差别，所以任意字段组合的结果都与完整解码一致。名字、质量值同样不需要参考；碱基与 MD / NM 需要。
- **对齐基准**：`refs/htslib`（启用 bzip2、lzma）在测试中按 `hts_opt_add` 的格式选项写出 3.0 / 3.1、九种块压缩方法、内嵌参考、无参考、多参考 slice、小 slice 和多 slice 容器的 CRAM，与 `sam_itr_queryi`、`sam_index_build3` 逐字段、逐字节比较。

**变异（BCF）**：`omicsweave-vcf` 的 `BcfReader`，以 htslib 的读法为准。

- **打开**：只读 BGZF 压缩的 BCF 与它的 CSI；VCF 文本不按区间读，先转换成 BCF（`ow convert --index` 同时写出 `.csi`）。索引由调用方给出，Python 层在旁边找 `<bcf>.csi`，缺失即报 `IoError`。
- **查询**：与 BAM 共用 `omicsweave-bgzf` 的 `BinIndex`（BAI / CSI 的载入与 `hts_itr_query` 的 chunk 计算）、`ChunkWalk`（`hts_itr_next` 的 chunk 遍历）和 `SeekReader`（按虚拟偏移读 BGZF 文本）；记录的 `[pos, pos + rlen)` 用存储的 `rlen` 判断重叠，与 `bcf_itr_queryi` 逐条一致。只解码重叠的记录。
- **基因型矩阵**：`genotypes(region, samples)` 一次给出区间内每条记录的位置、参考跨度（htslib 解码的 `rlen`）、ID、REF、逗号连接的 ALT（Arrow 字符串列）和所选样本的 `int8 [记录, 样本, 倍性]` 等位基因编号（`-1` 缺失、`-2` 无此等位基因），以及 `bool [记录, 样本]` 的 phase。倍性取区间内最大值；各记录先按自身倍性暂存，倍性都相同时直接作为结果，否则只重排一次。剂量、频率等由调用方从矩阵计算。
- **对齐基准**：htslib 的 `bcf_itr_queryi` / `bcf_itr_next` 与 `bcf_get_format_values`，覆盖 SNP、indel、多等位、带 `END` 的符号等位、混合倍性、缺失和无 GT 的记录，以及不同 `min_shift` 的 CSI。

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
| `CramToBam` | `.cram`（+ 参考） | `.bam`（+ `.bai` / `.csi`） | 未压缩流与 htslib 读 CRAM 写 BAM 逐字节一致（header 为 htslib 解析后重建的文本）；索引与 `samtools index` 同一内容 |
| `BamSort`（`sort`） | `.bam` | 按坐标排序的 `.bam`（+ `.bai` / `.csi`） | 与 `samtools sort --no-PG` 记录顺序、header、块划分一致，level 0 时逐字节一致；索引与 `samtools index` 同一内容 |
| `BedToBigBed` | `.bed`（bedGraph 亦可） | `.bigbed` + `.sidecar` | 与 `BigBedWriter` 逐字节一致 |
| `BedGraphToBigWig` | `.bedgraph` | `.bigwig` + `.sidecar` | 与 `BigWigWriter` 逐字节一致；libBigWig 读回的区间、逐碱基值、总汇总一致，zoom 汇总与区间一致 |
| `WigToBigWig` | `.wig`（fixedStep / variableStep / bedGraph 行） | `.bigwig` + `.sidecar` | 同上 |
| `CoolToHic` | `.cool` / `.mcool` | `.hic` + `.sidecar` | 与 `HicWriter::write_matrix` 逐字节一致 |
| `PairsToHic` | 4DN `.pairs`、Juicer short（纯文本 / gzip / BGZF） | `.hic` + `.sidecar` | 与 `HicBuilder` 逐条 `add` 同样的接触逐字节一致；straw 读回的计数一致 |
| `BedMethylToBigBed` | bedMethyl（modkit / ENCODE） | `.bigbed` + `.sidecar` | 全部字段原样保留，附该方言的 autoSql |
| `BedMethylToBigWig` | bedMethyl（modkit / ENCODE） | `.bigwig` + `.sidecar` | 区间与所选列的值一致；与 `BigWigWriter` 同一写出路径 |

**选择转换。** `Format` 列出 OmicsWeave 读写的全部格式（含 `.fai`、`.gzi`、`.csi`、`.bai`、`.sidecar`），同一格式的不同存储（plain / gzip / bgzip）各是一个值，名称与 Python 的 `File.format` 一致；`Format::detect` 按内容、`Format::from_extension` 按名称识别。`Conversion::between(from, to, src, dst)` 按格式对选择转换，源格式的任一存储都可以；`Conversion::infer` 两端都按扩展名判断；`Conversion::configure(&Options)` 应用选项，不适用于该转换的选项报错。`Report.files` 带每个文件的 `Format`；`Report.warnings` 列出转换跳过而不报错的内容（例如染色体列表外的行），Python 以 `omicsweave.convert.ConvertWarning` 发出，命令行打印到 stderr。

**分工。** Python（`omicsweave.convert`、`ow convert`）只做顶层组装：把路径或文件对象交给原生层，把返回的格式名组装成文件对象，命令行负责参数、进度与输出。识别格式（源按内容、目标按名称）、校验选项、预算缺省值、解析大小都在 Rust。

**代码结构。** `format`、`conversion`（选择、`plan` / `run`、运行上下文）之外，每种源格式一个模块（`fasta`、`vcf`、`sam`、`bamsort`、`cram`、`bed`、`wig`（bedGraph 与 wig）、`bedmethyl`、`cool`、`pairs`（4DN pairs 与 Juicer short）），共用部件在 `pipeline/`：`input`（压缩识别、并行解压、按行分批）、`output`（临时文件与原子改名、spill 目录）、`sort`（外排序）、`records`（VCF → BCF、SAM → BAM 共用的「文本行 → BGZF 记录 + 分箱索引」流水线）、`hic`（cool → hic 与 pairs → hic 共用的「按文件顺序的像素 → 矩阵、块、`.hic` 与 sidecar」）。bedMethyl → bigBed 复用 `bed` 的解析—排序—写出流程（换成 bedMethyl 的行解析与 autoSql），→ bigWig 复用 `wig` 的写出段。

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
| BAM → 排序 BAM | 解压、按批校验记录、分区排序、deflate | 切记录、header、归并、BGZF 分帧、建索引 |
| CRAM → BAM | 按 slice 解压各块、解码记录、编码 BAM 记录（参考由共享读取器加锁读取）、deflate | 读容器、定位 slice、BGZF 分帧、建索引 |
| BED → bigBed | 解析、分区排序、块压缩 | 染色体编号、zoom sweep、R-tree |
| bedGraph / wig → bigWig | 解析、分区排序、块压缩 | 解析 wig 段首（lead）、染色体编号、zoom 汇总、R-tree |
| cool → hic | chunk 读取与解码、像素定位、分区排序、块编码与压缩 | 矩阵记录、footer |
| pairs → hic | 解压（BGZF）、按行批解析与最细分辨率分箱、分区排序、按固定大小的像素块汇总到各较粗分辨率、块编码与压缩 | 读 header、切行、合并同一格的计数、矩阵记录、footer |
| bedMethyl → bigBed / bigWig | 解析与校验、按 mod / strand 选择、分区排序、块压缩 | 染色体编号、重叠检查、zoom 汇总、R-tree |

普通 gzip（非 BGZF）无法拆分，只能在单独的线程里流水线解压，不能并行。

**BGZF（`omicsweave-bgzf`）。** fasta、vcf、sam、bam 共用。`Writer` 按 htslib `bgzf_write` / `bgzf_flush_try` / `bgzf_flush` 的规则切块，所以未压缩块与 htslib 一致。块在执行器上压缩；写入时报告的位置是 `Logical` 偏移（块序号代替压缩地址），文件写完后由 `Blocks::resolve` 换算成虚拟偏移。`BinningIndex` 按 `hts_idx_push` / `hts_idx_finish` 的逻辑建 BAI / CSI：先在 `Logical` 偏移上做在线合并，偏移换算之后再执行那些需要比较地址大小的步骤（合并小 bin、合并同块 chunk）。`Reader` 在执行器上预读并解压后续块。读取端的 `BinIndex`（载入 BAI / CSI，按 `hts_itr_query` 算出区间的 chunk）、`ChunkWalk`（按 `hts_itr_next` 遍历 chunk）和 `SeekReader`（按虚拟偏移读文本，每个 chunk 一次 pread）由 BAM 与 BCF 的区间查询共用，记录格式各自解析。

**预算。** `Budget { threads, memory, tmp_dir }`。内存采用规划式分配，运行时不计量：`plan` 在读数据之前给出 `MemoryPlan`，包括三部分：
- `fixed`：常驻部分，如索引、header、bins 表；
- `pipeline`：每个队列的在途上限 × 单元大小；
- `buffer`：需要重排的转换使用的排序缓冲区。

`fixed + pipeline` 超出预算时直接报错，不在运行中途 OOM。批大小随预算自适应（64 KiB 到 4 MiB）。核算只覆盖大块缓冲，分配器开销、线程栈和小元数据不计入，RSS 允许超出约 10% 加 64 MiB。

**超过内存时。**
- **纯流式**（FASTA → bgzip、VCF → BCF、SAM → BAM）：工作集只有在途批次，任意大小都能处理。
- **单元级**（FASTA → 2bit）：2bit 记录头在打包数据之前，所以一条序列的打包数据先暂存；超出缓冲时写到临时文件。目录在所有记录之前，序列名来自源文件的 `.fai`，没有 `.fai` 时多扫一遍输入。
- **需要重排**（BED / bedMethyl → bigBed、bedGraph / wig / bedMethyl → bigWig、BAM 排序、cool / pairs → hic）：用通用外排序 `pipeline::sort::Sorter`。缓冲满时，并行排序各部分，归并写成 run 文件，最后与剩余的内存部分做 k 路归并。键唯一（BED、BAM 加输入序号，bedGraph / wig 加行号），所以结果与在哪里 spill 无关。BAM 排序至少需要 1 MiB 的排序缓冲区（`samtools sort` 每线程的下限），预算容不下时报错。cool → hic 按 (染色体对, 块号, y, x) 排序，之后用 hic 的流式接口（`begin_matrix` / `add_block` / `end_matrix`）逐块写出，整张矩阵不需要放进内存。pairs → hic 的键相同；同一格的多条接触就是一个像素的计数，`Spill::ABSORBS` 让排序器在写 run 文件时就把相等的键合并（`absorb`），最后的归并流里调用方再合并一次。计数是整数，用 `f64` 求和是精确的，合并发生在哪里都不影响结果；Juicer short 的 score 可以是小数，这时最细分辨率的键带上行号、较粗分辨率的键带上像素块号，排序器不会合并它们，每个格子的和总是按同一顺序求得，与线程数和 spill 无关。

spill 文件默认放在输出文件旁边（不用 `/tmp`）。所有输出先写到 `.partial`，整个转换成功后才改名；失败时清理。

**内存数据结构。** 每种格式只定义转换需要的流式单元，不定义整文件模型，也不对外提供操作接口：

| 数据模型 | 单元 |
|---|---|
| VCF / BCF | 复用 `Record` |
| SAM / BAM | 行到 BAM 记录的编码器（`BamEncoder`）；BAM 记录（`Record`） |
| FASTA / 2bit | `TwoBitPiece` / `TwoBitRecord` |
| Interval | `(chrom, start, end, rest)` |
| Contact | 按 `MatrixLayout` 定位后的 `Record` |

查询接口仍然只有"a 文件 → b 内存产物"一种；Python 层只暴露 `omicsweave.convert`（`run` / `plan` / `routes`），命令行是 `ow convert`。

**与参照工具的有意偏差。**
- VCF 中 header 未声明的名字报错。htslib 会动态补 header，但写出的 BCF 会损坏；需要补 contig 时开启 `add_missing_contigs`，先扫一遍 CHROM 列。
- SAM 中 RNAME 不在 @SQ 里时报错。htslib 只警告，并把读段当作未比对。
- BED 的染色体按首次出现顺序编号；bigBed 的 B+ 树按名字索引，读取不受影响。bedGraph / wig → bigWig 同样如此（`wigToBigWig` 按名字排序）。
- wig 的 fixedStep / variableStep 段写成 bedGraph 类型的数据块（`wigToBigWig` 保留段类型）；解码出的区间与值相同。bedGraph / wig 不要求按序（UCSC 要求已排序输入），`track` / `browser` 行在任意位置跳过，span 大于 step 造成的重叠报错；空区间、非有限值报错；不提供 `-clip`，越过染色体末端报错。
- wig 按批并行解析：一批开头尚未见到段首的行先记为 lead，由串行阶段按上一批结束时的段补全位置，每行只解码一次。
- cool 只取上三角（`square` 模式同样如此），不读 balance weight。
- `.hic` 不加 Juicer 的 `All` 染色体。
- pairs → hic 不做 Juicer `pre` 的任何过滤（MAPQ、pair type、同一酶切片段内的读段对）：每行算一次接触（short-with-score 算它的 score）。位置按原样分箱（与 `pre` 相同，1-based 的位置不减一），超出染色体长度的位置计在末端（与 `HicBuilder` 相同）。染色体列表外的行跳过，在 `Report.warnings` 里按染色体给出行数。只写 BP 分辨率；FRAG 需要酶切位点文件，不支持。缺省分辨率是 `pre` 的 9 个（2.5M 到 5k）；较粗分辨率由最细的像素汇总，必须是最细分辨率的整数倍，否则报错并指出是哪一个。每格计数先精确求和再存成 `f32`（`pre` 与 `HicBuilder` 直接用 `f32` 累加，计数超过 2^24 后才会不同）。分辨率按 `pre` 的习惯从粗到细列出、按此顺序写出，最细分辨率的像素先写到 spill 目录里的临时文件，最后写出。
- bedMethyl 按字段数分方言：18 列为 modkit `pileup`（第 10–18 列在旧版 modkit 中以空格分隔，同样接受），11 列为 ENCODE WGBS；混用报错。→ bigWig 时值原样保存（百分比保持 0–100，不从计数推导比例）；文件里有多个修饰代码（第 4 列）而没有给 `mod` 时报错并列出代码；所选条目仍有重叠时报错，并指出该加的选项（两条链重叠时提示 `strand`）。
- BAM 排序不加 `@PG` 行（相当于 `samtools sort --no-PG`）；`--index` 缺省写 `.bai`（参考序列放不下时写 `.csi`），`samtools sort --write-index` 缺省写 `.csi`。BAM 到 BAM 只有排序这一种转换，没有 `sort` 时报错；未排序的 BAM 不能只建索引。

**BAM 排序。** 排序键与 `samtools sort` 的坐标序（`bam_sort.c` 中的 `ks_radixsort` 和归并时的 `heap_lt`）相同：参考序列编号按无符号数比较（`tid` -1 即无参考序列的记录排最后），再比 `pos + 1`（取低 63 位），再比反链标志，其余按输入顺序。记录按原始字节搬运，不重新编码，只写回 htslib 读入时重算的 `bin`（有 CIGAR 的记录）。header 按 htslib 读入、解析、重建的规则改写：`sam_hdr_sanitise`（每行以 `@` 开头，补上最后的换行）、行格式校验、二进制参考序列列表中没有 `@SQ` 行的补在最后一个 `@SQ` 之后、`LN` 不一致时以二进制列表为准、只在文本中出现的 `@SQ` 追加到参考序列列表；第一个 `@HD` 移到最前，`SO` 改为 `coordinate`，删去 `GO`、`SS`，没有 `@HD` 时补 `@HD VN:1.6 SO:coordinate`。htslib 拒绝的 header（缺 `SN` / `LN` 的 `@SQ`、缺 `ID` 的 `@RG` / `@PG`、重复的参考序列等）同样报错。与 `samtools sort` 的差异只在病态输入：read name 缺结尾 NUL、CIGAR 存放在 `CG` 标签中的记录按原样复制（htslib 会补 NUL、把 `CG` 移到标签末尾或把短 CIGAR 写回原位）。参照实现 `refs/samtools`（samtools 1.25 的 `bam_sort.c` 及其依赖，进程内调用）只用于测试。

### 4.7 Python 层（`omicsweave`）

Python 绑定端到端处理加载：按路径打开文件、按 key 查询、拿到 numpy 结果，不暴露索引、数据块、解压器等中间数据。

- **类型都定义在 Python 里**：key（`Region`、`ContactWindow`、consensus 的 `Haplotype` / `Samples` / `AllAlt`）、结果（每个结果带上它的 key，数据全是 numpy 数组）、选项和异常都是 Python dataclass 或异常类；原生模块 `omicsweave._native` 只收发基本类型和 numpy 数组。
- **三层**：`types` 是与格式无关的 key 与结果；`formats` 每种文件格式一个类，文件对象只描述路径上的格式，同一格式的不同存储（plain / gzip / bgzip）各是一个类，索引与伴随文件（`.fai`、`.gzi`、`.sidecar`）也是平等的格式类，writer 属于各自的格式；`query` 是读取器，接收文件对象，按 key 返回 `types` 的结果。格式按内容判断（原生 `detect`，gzip 文件解压开头几个字节），与扩展名无关。
- **key 与方向**：`Region(chrom, start, end, strand=".")` 带链；`Spliced(chrom, spans, strand, name)` 是若干升序片段，两者合称 `Locus`，是稠密查询的 key（`SequenceReader.fetch`、`TrackReader.signal` / `bins`、`ContactReader.matrix` 的行列区间）。稠密结果带 `locus` 字段，稀疏结果带 `region` 字段。`Spliced` 提供 `to_genome` / `to_local`，在输出位置与基因组坐标之间向量化换算。`Annotation.locus(transcript, part)` / `loci(rows, part)` 由注释给出转录本的 `exons`、`cds`（含终止密码子）、`utr5`、`utr3`；`StrandedTrackReader(plus, minus)` 按 key 的链读一对分链 bigWig，无链 key 报错。原生层只接收 `(chrom, [(start, end)], strand)`。
- **列优先**：区间、track 区间、接触 cell 按列输出成多个连续数组，不转成 Python 对象；字符串列用 Arrow `LargeString` 布局（`offsets: int64`、`data: uint8`），可零拷贝转成 pyarrow。core 的行优先 entry 只在原生层内部暂存，拆列是那一次拷贝。
- **内存**：所有查询默认自动分配并直接返回结果。定长输出由 Python 或原生绑定层分配 numpy 数组，变长结果在原生层按查询得到的大小分配，交付 numpy；内部的描述、分配、读取、解码和写入步骤不拆成 Python 用户必须调用的接口。每次默认查询的结果独立持有数据，后续查询、缓存淘汰或关闭 reader 不改变已返回的数组。定长查询保留可选 `out=`，要求可写、C 连续、形状与 dtype 匹配，返回结果引用传入的同一数组；调用方再次写入该数组时，共享它的结果随之改变。
- **缓存**：普通查询无需创建 cache、描述依赖或预取。`cache=` 及其控制接口用于需要共享解码块或规划驻留的调用方；输出数组和解码工作区的自动管理不依赖缓存是否启用。
- **伴随文件只显式构建**：FASTA 的 `.fai`（bgzip 还有 `.gzi`）是按 key 查询的前提，可显式给出路径，否则在源文件同目录按同名查找；打开时从不自动构建，缺失即报错，由 FASTA 文件对象的 `build_index()` 手动构建或由 `FastaWriter` 写出。BAM 的 `.bai` / `.csi` 同理：按 htslib 的顺序在旁边查找，缺失即报错，由 `Bam.build_index()` 构建或在 SAM → BAM 转换时一并写出。按 key 直接查询不读 sidecar：从已有文件建 sidecar 要遍历整个文件，单次查询用不上；它只供以后的 loader 使用，由文件对象的 `build_sidecar()` 或 writer 的 `sidecar=True` 显式产出。
- **并发**：每个读取器内部状态由一把锁保护，只在释放 GIL 后加锁，解码期间不持有 GIL；读取器可 pickle，按文件和选项重新打开。
- **异常**：`Error` 的每个变体对应 `omicsweave.errors` 中的一个异常类，字段变成属性，并同时继承语义最近的内置异常。

### 4.8 注释（`omicsweave-gtf`）

`omicsweave-gtf` 把基因注释（GTF 2.2、GFF3）读成一个通用的、可复用的内存模型，供定量（BAM 读段按基因 / 外显子 / 内含子 / 剪接位点计数）、单细胞 ATAC 的 gene activity、按注释采样等功能共用。它不是训练数据加载的一部分：注释整体常驻（2.3 所说的常驻对象），按区间和按 id 查询。

**读取。**
- 支持 GENCODE、Ensembl、RefSeq、UCSC 导出的 GTF 和 GFF3；纯文本、gzip、BGZF 按内容识别。方言按 `##gff-version 3` 或第一条记录的属性语法识别，也可显式指定。
- 文本按行切批，在 `omicsweave_core::exec` 的线程池上并行解析，每批得到使用批内字符串编号的列；合并在单独的线程上按批次顺序进行，把批内字符串依次并入全局字符串池。所以字符串编号、行顺序与线程数和批大小无关，测试检查这一点。普通 gzip 在单独线程上解压，BGZF 在池上并行解压。
- 规范化：坐标统一为 0-based、左闭右开的 `u32`（只在读写时换算）；链为 `Strand`（`+`、`-`、`.`、`?`）；phase 为 0–2 或无；score 按原文保存（`.` 为空），因为 UCSC 写成 `0.000000`，转成浮点就不能原样写回。
- 层级是 gene → transcript → feature。GTF 按 (seqid, `gene_id`) 和 (seqid, `transcript_id`) 分组，同一个 id 出现在两条染色体上（PAR、UCSC 的多位置映射）就是两个基因；RefSeq 基因行的 `transcript_id ""` 视为没有。GFF3 按 `ID` / `Parent` 建树，转录本靠结构而不是类型名识别：有外显子级子记录（exon、CDS、UTR、密码子等）的记录是转录本（`mRNA`、`lnc_RNA`、`V_gene_segment` 都可以），转录本的父记录是基因；同一 `ID` 的多行外显子级记录（RefSeq 的 CDS）是一个不连续特征的多行；有多个父转录本的特征在每个转录本下各有一行。
- 缺失的 gene / transcript 记录由子记录合成：范围取子记录的最小 / 最大坐标，链取第一个子记录的链（UCSC 导出没有 gene / transcript 行；RefSeq 线粒体基因的 CDS 直接挂在 gene 下，合成一个转录本）。
- 名称、biotype、tag 从各来源的常用键读取：`gene_name` / `gene` / `Name`，`gene_type` / `gene_biotype` / `biotype`，`transcript_*` 同理；`tag` 可以重复（GTF）或逗号分隔（GFF3）。常见 tag（`basic`、`Ensembl_canonical`、`MANE_Select`、`CCDS`、`appris_principal_*`、`RefSeq Select` 等）同时存为位图，大小写和空格 / 下划线的差别不影响匹配，所有 tag 值另以字符串列表保存。这些字段与属性保留策略无关，总会读取。
- 染色体别名：读取时可按表改名；`match_seqids` 把另一个文件的名字（BAM header 的参考名）对应到注释的 seqid，依次试原名、别名表（双向）、内置的 UCSC ↔ Ensembl 写法（`chr1` ↔ `1`、`chrM` ↔ `MT`）。
- 校验默认只收集报告：孤儿记录、未知 `Parent`、重复 id、链冲突、seqid 冲突、子记录超出父记录、不合层级的嵌套、多个父转录本。`strict` 时第一个问题即报 `Error::Format`。输入不要求排序。语法错误（列数、坐标、链、phase、引号）总是报错。

**内存模型。** 不为每条记录分配堆对象：
- `StringPool`：一块字节 arena 加偏移，内容查找用开放寻址的 id 表。seqid、source、类型、id、名称、biotype、属性键和值都只存一次，用 `u32` 引用。
- `Genes`、`Transcripts`、`Features`、`Others` 是列式结构。一个基因的转录本、一个转录本的特征都是连续的行区间（CSR 偏移），按位置排序；基因按 (seqid, start, end) 排序，所以基因行号即基因组顺序。不属于层级的记录（GFF3 的 `region`、`chromosome`，有问题的记录）进 `Others`，保留文件顺序。
- 属性（第 9 列）按文件记录编号存成一张平表：每个不同的 (key, value) 对只存一次，记录只存对的编号。GENCODE 外显子行重复转录本的十几个属性，按对去重后，1.6M 行、全部属性的模型比直接存 (key, value) 小约 30%。保留策略 `All`（无损、原顺序，写回原 GTF 需要）、`Only(keys)`、`None`。GTF 值是否带引号记在键编号的最高位。
- 按 id、按名称的查找表首次使用时构建；索引、派生视图同样按需构建并缓存。注释构建完成后不可变，可在线程间共享。

**区间索引。** `IntervalIndex` 按 seqid 分组，组内按 (start, end, id) 排序存四个平行数组（start、end、id、子树最大 end），即 cgranges 的隐式增强区间树：树就是数组本身，构建是一次排序加一次自底向上的扫描，查询只访问连续内存。基因、转录本、外显子各有一个索引，id 是对应表的行号。查询：overlap（按 start 顺序给出）、contained_in、containing、nearest（重叠为 0，否则取两侧最近的全部并列项）、批量查询（CSR 输出）。结果追加到调用方的缓冲区。

**SweepCursor。** 坐标排序的 BAM 逐条扫描时，读段起点不减。游标维护「已进入」（起点小于某次查询终点）且「未过期」（终点大于下限 floor）的活动集合：每个项目只进入、过期各一次，推进是摊还 O(1)；一次查询的代价是活动集合的大小，不分配内存。剪接读段的多个比对块彼此不单调，所以分成两步：每条读段 `advance(read_start)` 一次，再对每个块调用 `overlaps`；起点不低于 floor 的任何查询都精确。低于 floor 的查询从树重新定位（O(log n + k)），结果仍然正确。

**派生视图。** 按需计算，无参数的视图缓存在注释上：
- 每个转录本的外显子块（相交或相接的外显子合并；没有 exon 行时用 CDS / UTR / 密码子行），内含子，基因的合并外显子区，注释剪接位点集合（每个不同的 (seqid, start, end, strand) 一次，附 donor / acceptor、所属基因和转录本，可按坐标精确查找），DEXSeq 式不相交外显子 bin（按基因，或按同链外显子重叠的基因聚合），转录本长度、CDS 长度。
- 计数特征窗口：TSS / TES，promoter（`GenomicRanges::promoters` 的定义）和 TES 窗口（上下游可配置、按链方向），基因体加上游延伸（ArchR / Signac 的 gene activity）；在 0 处截断，给出染色体长度时在末端截断。没有链信息的记录按 `+` 处理。
- 区间视图都是 `IntervalSet`：(seqid, start, end, strand, owner, name) 列加自己的索引和游标，owner 是所属基因或转录本的行号，所以计数层对它们和对基因、外显子用完全相同的方式索引和扫描。peak、固定 bin 等其它来源的特征也用这个类型（`IntervalSetBuilder`、`IntervalSet::from_columns`）。
- `Derive` 在 `Selection` 上、按染色体长度截断地构建这些视图。`select` 按 biotype、tag、基因 id / 名称、染色体给出基因和转录本的行号列表，不复制数据。

**输出。** GTF 有两种布局：`Original` 按文件顺序写回全部记录和 `#` 行，对全部属性读入的 GTF 逐字节一致（合成的记录不写）；`Model` 按模型写 gene、transcript、feature 行，先写模型字段（`gene_id`、`transcript_id`、名称、biotype、tag），再写其余保留的属性，GFF3 也由此转成 GTF。BED12 按转录本写，block 是外显子块，thick 部分覆盖 CDS 和密码子，非编码转录本按 UCSC 惯例 `thickStart = thickEnd = chromEnd`。

**Python。** `omicsweave.annotation.Annotation.read(path, attributes="all" | [keys] | None, aliases=None, strict=False)`：`genes`、`transcripts`、`features`、`exons` 是列式 dataclass（numpy 列，字符串为 `StringVector`），`gene(id)`、`overlaps(region, level)`、`nearest`、`subset(...)`、`introns()`、`junctions()`、`exon_bins()`、`promoters()`、`gene_bodies()`、`match_chroms()`、`report`。原生层只收发基本类型和 numpy 数组；注释不可变，不需要锁，调用时释放 GIL。pickle 时按路径和选项重新读取。

**性能。** 合成的 GENCODE 规模文件（1.6M 行、458 MiB、47k 基因、207k 转录本，`examples/parse_bench.rs`）：单线程读取 3.4 s，4–8 线程 0.9–1.3 s（随机器负载变化），全部属性的模型 228 MiB；不保留属性时 0.66 s、108 MiB。三个索引共约 0.07 s，全部派生视图约 0.45 s；外显子游标每次块查询约 7 ns。

**转录本的部分。** `Annotation::part_blocks(t, TranscriptPart, out)` 给出转录本某一部分的升序合并片段，作为 locus 的 spans：

- `Exons`：外显子块，即成熟转录本。
- `Cds`：CDS 与密码子行合并，终止密码子包括在内。GENCODE / Ensembl GTF 把 `stop_codon` 列在 `CDS` 之外，GFF3 通常把它算进 CDS，合并后两种来源结果一致，完整 ORF 的长度是 3 的倍数。
- `Utr5` / `Utr3`：外显子块中位于编码块 5' / 3' 一侧的部分，由外显子减去编码区间得到，不依赖 `UTR` 行有没有区分 5' 和 3'（GENCODE GTF 不区分）。非编码转录本没有这两部分。

拼接后的序列和信号由 Sequence / Track 读取器按 locus 读出（见 4.3），注释层不提取序列。

**不做。** 不做 gffread 式的转录本合并、去冗余；不解码 GFF3 的百分号转义（值按原文保存）；不写 GFF3。

### 4.9 定量层（`omicsweave-quant`）

`omicsweave-quant` 把已比对、按坐标排序的 BAM 一次扫描成覆盖度轨道、计数矩阵（含 velocity 与标记层）、单细胞 ATAC fragments、剪接位点表、代谢标记统计、DNA 甲基化、Hi-C 接触和 SNV 等位基因计数。边界：不做比对，不链接、不复用 STAR；与格式 crate 和 `ow convert` 不同，这一层做分析（过滤、特征分配、计数、去重），1.1 中「转换不做业务逻辑」的约束不适用于它。

**一次扫描，多个状态机。** 每个产物是一个累加器（`Accumulator`），所有累加器共享同一次扫描：

1. `source`：打开时检查 header 的 `@HD SO:coordinate`，扫描时检查 (tid, pos) 单调不减；任一不满足即报错，并提示 `ow convert in.bam sorted.bam --sort`。BGZF 块在线程池上解压（`omicsweave_bam::Scanner`），读取线程把连续记录切成批（与 convert 相同，64 KiB 到 4 MiB，随预算变化），遇到第一条未放置的记录（tid = -1，按序排在最后）即结束。多个 BAM 依次扫描，各自是一个独立的行来源。BAI / CSI 只会用于限定区域，目前不读。
2. `prelude`：`ReadView` 在第一次被问到时才计算并缓存派生量，同一条记录的多个累加器共享：FLAG / MAPQ / `NH` 过滤（`Filter`）、可选字段查找、文库链特异性、比对块（`M`/`=`/`X`，可选并入 `D`）、内含子（`N`）、5′ / 3′ 端、片段（由最左 mate 的 TLEN 得出）、Tn5 校正片段和切点、比对段（`segments()`：每个 `M`/`=`/`X` 的参考起点、读段起点、长度，供 `query_at` / `reference_at` 双向查位置）、读段下的参考碱基（`reference_bases()`：比对区间两侧各多取 2 bp，供序列上下文）、读段区间内的已知变异位置（`masked()`），以及对齐碱基（`aligned()`：每个 `M`/`=`/`X` 位置的参考位置、读段位置、读段碱基、参考碱基、质量、是否落在已知变异上）。参考碱基来自任务的参考基因组（`reference`：bgzip / 普通 FASTA 加 `.fai` / `.gzi`，或 2bit，经 `omicsweave_fasta`），没有参考或该染色体不在参考里时用 `MD` 标签；已知变异来自 BCF + CSI（`omicsweave_vcf::BcfReader`）。参考和变异在计划时按名字（含别名与 `chr1` ↔ `1`）对到每个输入的染色体，长度不符即报错；每批一个只向前滑动的窗口（按记录位置取序列、查询变异，前后预取 256 kb），每个碱基每批约取一次。mate 配对分两级：批内按 (read name, 本端位置, mate 位置) 配对（`Pairs`），配不上的半对在合并时交给 `Mates`；当合并前沿越过 mate 应在的位置仍未出现时，该半对按单端计。
3. `assign`：特征空间（`FeatureSpace`）是一个 `omicsweave_gtf::IntervalSet` 加上「行 → 列」映射：基因的合并外显子、基因全长、启动子 / gene activity 窗口（`Derive` 按 BAM 的染色体长度截断）、bigBed 的全部 peak（BED 不直接读，先用 `ow convert` 转成 bigBed）、固定大小的 bin。任务内按构造方式去重注册；每批每个空间一个 `SweepCursor`，染色体名用 `match_seqids` / `match_names` 对应（`chr1` ↔ `1`、别名表）。一条记录对同一空间、同一查询（块、片段、切点等）、同一链要求只分配一次。
4. `accum`：见下文。
5. `model` 与 `sink`：产物是与文件格式无关的数据（`Track`、`CountMatrix`、`FragmentTable`、`BarcodeMetrics`、`JunctionTable`、`Histogram`、`MismatchTable`、`Cytosine` 流、`Mbias`、`ContactTable`、`VariantTable`；矩阵的列是 `FeatureList`，特征表或按需生成行的位点列表），大的产物是只读一次的有序流（记录可能在 spill 文件里），由 sink 写成文件。

**Accumulator。** 每个累加器有两种状态：批状态（`Batch`）在 worker 上对一批连续记录 `observe`，结束时 `seal`（排序、折叠、批内配对，仍在 worker 上）；总状态（`Total<'p>`）在单个合并线程上按批次顺序 `merge`，合并时得到该批最后一条记录的坐标键（frontier），据此判断哪些结果已经最终确定。`finish` 把总状态变成产物，`deliver` 交给 sink 并返回 `summary.json` 中的一段。类型擦除的 `Stage` / `Merging` 让任务持有任意组合。与最初草案的差别：没有单独的 `spill(state)` / `finish(state, runs)`，批状态和总状态分成两个类型，总状态在 `RunCtx` 给出的预算份额内自己 spill（外排序 `Sorter`、有序缓冲 `Spool`）。`Needs`（cigar、tags、mates、assign、sequence、quality、reference、variants）声明累加器读什么，计划时取并集；只有并集含 reference / variants 时批才建参考 / 变异窗口，计划时为每个在途批多计 2 MiB / 1 MiB。

新增一个产物只需写一个累加器：声明 `Needs`、`estimate` 内存（不可 spill 的固定部分、每批在途部分、最小缓冲和缓冲权重）、批如何观察一条记录、批如何合并、交付什么。

**预算与 spill。** 复用 convert 的 `Budget { threads, memory, tmp_dir }` 和 `MemoryPlan`：`plan` 只读 BAM header 并构建特征空间，给出 `fixed`（特征表、索引、白名单、注释）、`pipeline`（在途的块和批）、`buffer`（各累加器的排序缓冲和 spool）。`fixed + pipeline + 各累加器最小缓冲` 超出预算时直接报错；其余缓冲按权重分给各累加器。总状态超出份额时 spill 到 `tmp_dir`（默认输出目录），不能 spill 的部分（无白名单时的条形码字典、等待配对的 mate）超出份额时报错并提示给白名单或增加内存。所有输出先写到临时文件，全部产物成功后才改名；产物在各自线程上并行写出，一起发布。

**确定性。** 批次按序合并，所有折叠满足交换律与结合律（整数求和、取最大、去重计数）；行号按条形码首次出现的记录顺序编号（与批的切分无关）；覆盖度的 run 跨批合并成极大段。所以输出与线程数、批大小、是否 spill 都无关，测试逐字节检查（`summary.json` 只记录与此无关的统计量）。

**覆盖度（`Coverage`）。** 单位：比对碱基（`samtools depth` 的语义，可选计入 `D`，即 `-J`；`N` 不计）、整段 span、5′ 端、3′ 端、片段（每对只从最左 mate 计一次，`[pos, pos + |TLEN|)`，可设最大长度；单端读段按 span，非最左 mate 不计）。链：合并，或按文库类型拆成正负两条（read 2 视作 read 1 的链，反向文库再翻转）。过滤默认排除 `0x704`，可设最小 MAPQ 和只计 `NH == 1`。缩放：原值、CPM（每百万计入的读段或片段）或给定因子。批把区间端点记成排序的事件；总状态保留前沿之后的事件，把之前的扫成等值 run 写入 spool，内存只随前沿处的深度增长，与基因组大小无关；最后写成 bigWig（`BigWigWriter`，带 zoom，可选 sidecar）。

**计数矩阵（`Counts`）。** 一个累加器，维度可组合：

- 特征（列）：上面的任一特征空间。
- 单位：读段；read pair（两端比对块的并，featureCounts `-p --countReadPairs`；不带 `multi_overlap` 时按 featureCounts 的投票：被两端都覆盖的特征优先，并列才算模糊）；片段（最左 mate 的 TLEN 区间，可做 Tn5 校正）；Tn5 切点（校正片段的首、末碱基，各自分配）；分子 = 读段 + UMI 去重。
- 行：标签（`CB` 或任意标签）、read group（各 header 的 `@RG` 依次合并）、输入文件；可给白名单（行按白名单顺序，含没有读段的条形码）。
- 去重：无；UMI（每个 (行, 特征) 中不同的 `UB` 精确值计一次，跨批全局去重：条目带 UMI 进入外排序，读出时折叠）；位置（同一行、同一片段位置只计一次，ATAC 的重复片段）。
- 分配规则同 featureCounts：至少 1 bp 重叠；链特异性 0 / 1 / 2；命中多个列的单位记为模糊（`multi_overlap` 时计入每一列）；多重比对（`NH > 1`）默认不计。

批输出单元格条目 `(row, col, layer, sub, value)`（`sub` 是 UMI 或片段位置），在批内排序折叠；总状态把条目送进外排序，矩阵按行读出时再按 `Fold` 折叠：计数层和标记层按 `Reduce`（求和、按 sub 去重计数、按 sub 取最大），velocity 证据层按位或后分类到三层；内存由排序缓冲界定，与非零元数量无关。

预设：单细胞 RNA（genes × UMI × `CB`，STARsolo `Gene` 的原始矩阵：全部条形码、不做细胞判定、`CB` / `UB` 按原样使用、正链文库、只计唯一比对、排除 `0x704`）；bulk RNA（genes × read 或 pair × 文件，featureCounts 默认：除未比对外全部记录都计、多重比对按 `NH` 排除，另写 featureCounts 格式的 `counts.tsv`）；scATAC（peaks / tiles × 切点 × `CB`，gene activity（基因体加上游 2 kb，Signac）× 片段 × `CB`，均按位置去重）；ChIP / CUT&Tag bulk（bins 或 peaks × 片段 × 文件，片段计入它重叠的每个特征）。

**单细胞 ATAC（`Atac`）。** 定义（10x Genomics `fragments.tsv`、ArchR / Signac）：

- 片段：通过过滤（默认要求 `0x1 | 0x2`，排除未比对、mate 未比对、secondary、QC fail、supplementary，MAPQ ≥ 30；重复保留）、两端比对到同一染色体的 read pair，取最左 mate：`[pos + 4, pos + |TLEN| - 5)`，即最左 mate 5′ 端 +4 到另一 mate 5′ 端 -5。
- 切点：片段的首碱基和末碱基（`start`、`end - 1`）。
- 重复：同一条形码、同一 (染色体, start, end) 是一个片段，第五列是其 read pair 数。
- 每条形码指标（`barcode_metrics.csv`）：read pair 数、片段数、有切点落在 peak 内的片段数、peak 内切点数、TSS ±50 bp（101 bp）内切点数、距 TSS 1900-2000 bp 的两个 100 bp 侧翼内切点数；ArchR 的 TSS enrichment 为 `(tss_cuts / 101) / (tss_flank_cuts / 200)`。

片段在批内排序折叠，总状态把前沿之前的片段（位置不可能再有重复的）按序写入 spool，最后写成 BGZF 的 `fragments.tsv.gz`，同时用 `BinningIndex`（BAI 参数）建索引，换成 tabix 的头部（`tbx_conf_bed`、按出现顺序的序列名）得到 `.tbi`。ATAC 的各矩阵与片段文件使用相同的过滤和片段定义，因此各目录的 `barcodes.tsv.gz` 一致。

**剪接（`Junctions`）。** 定义（STAR `SJ.out.tab`）：

- 剪接位点是 `N` 操作的参考区间。过滤默认排除未比对、QC fail、supplementary（STAR 的嵌合片段），重复和 secondary 保留；多重比对的每条比对都计。
- 一对 read 同时跨过同一位点只计一次（`pairs`，按 read pair 计数同样的配对；配不上的半对单独计）。
- unique / multi：`NH` 为 1 / 大于 1 的跨越读段（对）数。
- overhang：读段在该位点两侧相邻比对段中较短的一段，比对段是不被插入、缺失、clip、内含子打断的 `M`/`=`/`X`；表中取最大值。
- motif：参考序列上内含子首两个与末两个碱基：1 `GT/AG`、2 `CT/AC`、3 `GC/AG`、4 `CT/GC`、5 `AT/AC`、6 `GT/AT`，其余（及无参考时）为 0；链：奇数 motif 为 `+`，偶数为 `-`，非经典位点有注释时取注释的链，否则 0。annotated：等于注释中某个转录本的内含子（`Annotation::junctions()`，`find_any`）。

批输出 `(tid, start, end, unique, multi, overhang)` 条目，批内排序折叠，总状态进外排序；结束时折叠成 `JunctionTable`（按位置排序，motif 在合并线程上逐位点取 4 个碱基），写成 `junctions/junctions.tab`（STAR 的九列）。可选的 junction × 行矩阵（`JunctionMatrix`：行键、白名单、UMI 去重、额外过滤；STARsolo `SJ` 预设为 `CB` × UMI、只计唯一比对）：列是表中全部位点（按表顺序，含只有多重比对读段的位点），条目先按位点排序（列号在所有位点已知后才确定），读出时编号再进第二个外排序按行读出；与表写在同一个 `junctions/` 目录下。多个 BAM 要求参考序列相同且顺序一致。

**Velocity（`Counts::velocity`）。** 基因矩阵上的 spliced / unspliced / ambiguous 三层（velocyto / STARsolo `Velocyto` 规则）：读段的比对块与恰好一个基因的全长（含内含子，按文库链）重叠才归属该基因，重叠多个基因不计；对该基因每个包含读段全部比对块的转录本，每个块都落在同一外显子内为 exonic，否则（块在内含子里或跨外显子-内含子边界）为 intronic；没有转录本包含读段则不计。读段的证据位（exonic、intronic）按分子（有 UMI 时 `(行, 基因, UMI)`；无 UMI 时单个读段，read pair 按 read name 哈希）在 `value` 中按位或（`Fold`：该层按位或，其余层用 `Reduce`），读出时只有 exonic 为 spliced、只有 intronic 为 unspliced、两者都有为 ambiguous。与计数层的基因分配相互独立（计数层仍按外显子）。

**代谢标记（`Labeling`、`Counts::labeled`）。** SLAM-seq / scSLAM-seq / TimeLapse / 4sU：

- 读段的转录链按文库：read 2 视作 read 1 的链，反向文库翻转（`Labeling` 把非链特异文库按正链读；基因矩阵在非链特异文库下用所计基因的链）。转换按转录本方向给出（默认 `T>C`，可配置任意 `X>Y`）：正链转录本在基因组上是 `T>C`，负链是 `A>G`。
- 只计对齐碱基中参考与读段碱基都是 `A`/`C`/`G`/`T`、质量不低于 `min_quality`（默认 20；没有质量的记录全部通过）、不在已知变异上的位置；覆盖度是转录本方向为源碱基（`T`）的位置数，转换数是其中读成目标碱基（`C`）的位置数。`convert` 一次算出两条链的结果（转换位置列表），供各累加器共用。
- 基因矩阵的 `labeled` 层：计数层计入的单位（读段、read pair、UMI 分子）中转换数不少于 `min_conversions` 的；read pair 合并两端的转换位置（同一位置计一次），分子只要有一条读段达到阈值就算标记。条目与计数层同 `sub`、同 `Reduce`，在单独的层号上，所以 `matrix.*` 就是 total。
- `Labeling` 产物：每读段 `(T 覆盖度, T>C 数)` 直方图 `labeling/conversions.tsv.gz`（稀疏三列）；`labeling/mismatches.tsv`：按 mate、测序 cycle（从读段测序 5′ 端起，含 hard clip）、转录本方向参考碱基，统计四种读段碱基；`summary.json` 给出总转换率和十二种替换的总体比率。没有参考也没有 `MD` 的读段计入 `unassigned_no_reference`。

**DNA 甲基化（`Methylation`）。** WGBS、EM-seq、单细胞亚硫酸氢盐、带碱基修饰标签的长读段（Bismark / MethylDackel / modkit 的定义）：

- 调用（call）是读段上的一个胞嘧啶：参考位置（`C`，反链胞嘧啶取 `G` 的位置）、链、上下文、状态。来源按记录选择（`CallSource::Auto`）：有 `XM` 用 Bismark 的 `XM`（`Z/z` CpG、`X/x` CHG、`H/h` CHH，大写为甲基化，`U/u` 与 `.` 不是调用；链取 `XG`：`CT` 正链，`GA` 反链）；否则有 `MM`/`ML` 用碱基修饰（`C+m` 列表里每个 C 的概率 `(ML + 0.5) / 256` 不小于 `min_probability` 为甲基化、不大于 `1 - min_probability` 为未甲基化，中间不调用；未列出的 C 在默认 / `.` 模式下为未甲基化，`?` 模式下不调用；`MN` 与读段长度不符或有 hard clip 而无 `MN` 时视为标签不合；链取读段的）；否则对参考比对（bwa-meth / MethylDackel）：原始正链读段（`YD:Z:f`、`XG:Z:CT`，都没有时 read 1 正向 / read 2 反向）在参考 `C` 上读到 `C` 为甲基化、`T` 为未甲基化，原始负链读段在参考 `G` 上读到 `G` / `A`，用 `aligned()`。
- 上下文取参考（有参考基因组用基因组，否则用 `MD` 推出的、同一读段相邻对齐位置的参考碱基；`MM` 调用没有参考时用读段自身碱基）：其后为 `G` 是 CpG，否则 `H` 后接 `G` 是 CHG，否则 `HH` 是 CHH（反链在互补链上向左看）；后两位未知的胞嘧啶不调用。只调用质量不低于 `min_quality`（默认 5；没有质量的记录全部通过）、不在任务已知变异上的碱基。
- 重叠 mate 只计一次（Bismark `--no_overlap`）：同一染色体、相距小于 64 kb 的一对在 `Pairs` / `Mates` 中配对（两端判断相同，配不上的半对过期后按单端计），read 2 落在 read 1 比对区间内的调用丢弃。
- M-bias：按上下文、mate、测序方向的读段位置（1-based，含 hard clip）统计甲基化 / 未甲基化调用，在丢弃重叠之前；写 `mbias.tsv`，并写入 `summary.json`。
- 每胞嘧啶表与轨道：批内按 (tid, pos, code) 折叠，总状态保留前沿之后的胞嘧啶，并且不越过等待配对的半对中最小的调用位置（这些半对的第一个调用键放在一个计数多重集中），之前的写入 spool；最后一次读出同时写 `cytosines.tsv.gz`（chrom、start、end、context、strand、methylated、unmethylated，0-based 半开，BGZF，tabix 按 BED 建索引，与 ATAC fragments 共用 `IndexedBed`）、`beta.bw` 和 `coverage.bw`；`merge_cpg` 把一个 CpG 两条链合成一行 `[C, C + 2)`、链为 `.`。只支持单个 BAM（关掉表和轨道后矩阵可以跨多个 BAM，如每细胞一个 BAM）。
- 区域矩阵（`MethylMatrix`，allcools / scBS-map 风格）：特征为任一特征空间（基因、启动子、gene activity、peak、bin），行键同计数矩阵；每个调用按位置分配（`assign_at`，同一批的扫描游标），条目 `(row, col, layer, 0, 1)` 以 `Fold::direct(Sum)` 折叠成 `methylated` / `total` 两层（各层文件都按层命名：`methylated.cpk`、`total.cpk`）。

**Hi-C / Micro-C（`Contacts`）。** 输入是按坐标排序、带 mate 标签（`MC`、`MQ`，`samtools fixmate -m`）的 BAM（4DN pairs、pairtools 的定义）：

- 接触取自一对中坐标序在前的 mate（同位置取 read 1），用它的 RNEXT/PNEXT、mate 链（0x20）、`MC`（反向 mate 的 5′ 端需要它的 CIGAR）和 `MQ`，不缓存 mate。两端是两个 mate 的 5′ 位置（1-based），按染色体序号再按位置排成 `(chr1, pos1) <= (chr2, pos2)`。
- 过滤：默认要求 paired，排除未比对、mate 未比对、secondary、QC fail、supplementary、duplicate；两端 MAPQ 都不低于 `min_mapq`（默认 30，mate 的取 `MQ`，没有 `MQ` 只在阈值为 0 时通过）；同一染色体上距离小于 `min_distance` 的丢弃；`self_ligation` 时距离小于它的 inward（`+ -`，dangling end）和 outward（`- +`，self-circle）接触丢弃，同向接触不可能是自连接产物，保留。
- 去重：两端位置、链（和条形码）相同的接触只留文件序第一个。
- 批输出变长的 `Contact`（含 read name）进外排序，键 (chr1, chr2, pos1, pos2, strand1, strand2, 行, 输入, 记录号)，所以重复相邻且顺序与批切分无关；读出时去重，同时写 `contacts.pairs.gz`（BGZF，头部含 `#chromsize`，按 BAM 参考序）、把每个接触按每个分辨率分 bin（与 Juicer `pre` / `HicBuilder` 相同：位置截到染色体长度、`bin = pos / res`、染色体内取上三角）送进各分辨率自己的外排序（键 (c1, c2, block, y, x)，`ABSORBS` 折叠计数），最后粗到细逐张 map、按 block 序经 `HicWriter` 的 `begin_matrix` / `add_block` / `end_matrix` 写出 `contacts.hic`（v9，原始计数，可选 sidecar）。convert 的 pixel → `.hic` 阶段不是公开 API，所以这里直接调用 hic writer，业务逻辑（过滤、去重、两端定义）都在 quant。
- 单细胞：按 `barcode_tag` 统计每个条形码的 pairs、duplicates、contacts、cis、trans、cis ≥ 1 / 10 / 20 kb（pairtools stats），写 `cell_metrics.csv`，同时写入 `summary.json` 的 `cells`。

**等位基因计数（`Alleles`）。** cellsnp-lite mode 1 / vartrix：

- 位点是 BCF（带 CSI）中的二等位 SNV（REF、ALT 都是单个 `ACGT`），其他记录跳过；计划时整体读入（每个位点 12 字节，`SitePanel`），按文件序全部作为列（列表行按需生成，不为每个位点存字符串：`SiteFeatures` 实现 `FeatureList`）。
- 读段在位点上有质量不低于 `min_quality`（默认 20）的对齐碱基即观测到它：等于 REF 为 `ref`，等于 ALT 为 `alt`，其他 `ACGT` 为 `other`，`N` 不计；位点落在缺失或内含子中不计。
- 单位：UMI（`umi_tag` 精确值，没有 UMI 的读段不计）或 read name（一对的两端只计一次），即 `(行, 位点, 等位, sub)` 上的 `Reduce::Distinct`；一个单位的读段不一致时在它出现的每个等位上各计一次。行键同计数矩阵。
- 每位点表 `variants.tsv.gz`：矩阵写出时按列累加各层（`write_with` 的 `column_sums`），给出 REF / ALT / other 单位数和 `alt / (ref + alt)`。
- 基因水平 ASE（`AlleleGenes`）：取一个样本（默认第一个）定相杂合（`0|1`、`1|0`）位点，读段读到单倍型 1 的等位为第一位证据、另一等位为第二位；读段的比对块与恰好一个基因的外显子重叠（按文库链）时，证据按单位在 `(行, 基因, sub)` 上按位或，读出时分成 `hap1`、`hap2`、`both` 三层（与 velocity 共用两位证据的折叠）。

**输出布局。** 输出目录下：`<name>.bw`（或 `.forward.bw` / `.reverse.bw`）；计数产物 `<name>/`：`barcodes.tsv.gz`、`features.tsv.gz`（id、name、type、chrom、start、end 六列，Cell Ranger ARC 的列）、`matrix.cpk`（CellPack，`u32`，每行一条记录，按行顺序追加；默认）或 `matrix.mtx.gz`（Matrix Market，特征 × 条形码，同 Cell Ranger）或 `matrix.csv.gz`（稠密，一行一个特征），其他层为 `<layer>.cpk` 等（`spliced`、`unspliced`、`ambiguous`、`labeled`），bulk 另有 `counts.tsv`；`atac/fragments.tsv.gz`（+ `.tbi`）、`atac/barcode_metrics.csv`、`atac/peaks/`、`atac/tiles/`、`atac/gene_activity/`；`junctions/junctions.tab`（有矩阵时同目录为 10x 目录，特征类型 `Splice Junction`，id 为 `chrom:start-end`，0-based 半开）；`labeling/conversions.tsv.gz`、`labeling/mismatches.tsv`；`methylation/cytosines.tsv.gz`（+ `.tbi`）、`methylation/beta.bw`、`methylation/coverage.bw`、`methylation/mbias.tsv`、`methylation/<矩阵>/`（`methylated.*`、`total.*`）；`hic/contacts.pairs.gz`、`hic/contacts.hic`、`hic/cell_metrics.csv`；`alleles/`（`ref.*`、`alt.*`、`other.*`）、`alleles/variants.tsv.gz`、`alleles/genes/`（`hap1.*`、`hap2.*`、`both.*`）；`summary.json`（输入、`Needs`、每个产物的计数与各类未分配原因、spill 字节数、文件列表）。压缩文本都是 BGZF（与 gzip 兼容），在线程池上压缩。CellPack 通过 git 依赖（固定 rev）接入，不修改其代码。

**对齐测试。**

- 覆盖度：`refs/samtools` 中 vendor 的 `bam2depth.c`（`samtools depth -aa`，含 `-J`、`-Q`），逐碱基比较；其余单位、链拆分对照测试中按定义直接计算的逐碱基计数；1 线程与 8 线程、spill 与不 spill 的 bigWig 逐字节一致。
- bulk 计数：`refs/subread`（Subread 2.2.1 的 featureCounts，GPL-3.0，`publish = false`，只作测试依赖，不进入 wheel）在随机合成的 BAM 和 GTF（重叠基因、两条链、剪接读段、多重比对、重复、secondary / supplementary）上逐基因比较，单端 `-s 0/1/2`、`-O`，以及 `-p --countReadPairs`；`Assigned`、`Unassigned_NoFeatures`、`Unassigned_Ambiguity`、`Unassigned_MultiMapping`、`Unassigned_Unmapped` 也一致。
- 单细胞、read group、白名单、bins、三种矩阵格式：对照按定义直接计算的结果。
- ATAC：片段文件逐行、peak / tile 切点矩阵、gene activity、指标对照按上述定义的直接计算；`.tbi` 与 htslib `tbx_index_build` 对同一文件建出的索引结构一致（bin 顺序除外：htslib 按哈希表顺序写）。
- 剪接：随机合成的读段（成对、多重比对、重复、带插入缺失 clip 的剪接读段，另加正好剪在注释内含子上、两端都跨同一位点的 read pair），基因组随机生成后在大部分位点上写入六种经典 motif；`junctions.tab` 逐行对照按上述定义的直接计算（参考由 htslib `bgzf` + `fai_build3` 写成 bgzip FASTA，另写 2bit，两者输出逐字节一致；无参考时对照无 motif 的版本），junction × `CB` 矩阵和条形码顺序、velocity 三层（单细胞 UMI 与 bulk 读段）同样对照直接计算；1 线程与 7 线程、spill 与不 spill 的文件逐字节一致。
- 代谢标记：从随机基因组按 CIGAR（插入、缺失、内含子、soft / hard clip）抽取读段，正负链转录本、单端与 read pair，注入 `T>C`（转录本方向）、测序错误、`N`、样本 SNP 和随机质量；直方图、mismatch 表、单细胞（UMI，阈值 2）与 bulk（非链特异，按基因链）的 `labeled` 层对照直接计算，带与不带 SNP 掩码（BCF 由 htslib 写出并用 `bcf_index_build3` 建 CSI）；bgzip FASTA、2bit、`MD` 标签（测试按基因组生成）三种参考来源及不同线程数的输出逐字节一致。

- 甲基化：从生成的基因组抽取定向文库的 read pair（两条原始链、每个分子逐个胞嘧啶随机甲基化、测序错误、`N`、插入缺失与 clip、重叠 mate、被过滤的记录、已知 SNP），参考调用、`MD` 调用、`XM` 调用的 `cytosines.tsv.gz`（含 CpG 合并）、beta / coverage bigWig、M-bias、bin × `CB` 与 bin × 文件矩阵对照按上述定义的直接计算，`.tbi` 与 htslib `tbx_index_build` 一致；`MM`/`ML` 长读段（5hmC 列表在前、显式 / 隐式模式、两个概率阈值、有无参考）同样对照直接计算；1 线程与 6 线程、spill 与不 spill 的文件逐字节一致。
- Hi-C：生成的 read pair（cis / trans、两条链、clip 与插入缺失、`MC` / `MQ` 偶尔缺失、任一端低 MAPQ、带标志的记录、未标记的重复、条形码）；`contacts.pairs.gz` 与每条形码统计对照直接计算，`contacts.hic` 与 `HicBuilder` 对同一组接触写出的文件逐字节一致，并用 `HicReader` 读回 25 kb 的每张 map 对照直接分 bin；距离过滤、不去重、MAPQ 0 的 bulk 设置同样；不同线程数、spill 的文件逐字节一致。
- 等位基因：生成的 read pair（定相样本两条单倍型的分子、PCR 拷贝、错误、`N`、低质量、位点上的内含子与缺失、被过滤的记录），BCF 由 htslib 写出（含要跳过的插入和多等位记录）；单细胞 UMI 与 bulk read pair 的 ref / alt / other 层、`variants.tsv.gz`、hap1 / hap2 / both 基因层对照直接计算；不同线程数、spill 的文件逐字节一致。

**与参照工具的有意偏差。**

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

**性能。** 合成的单细胞文库（24 条 50 Mb 染色体，1200 万 read pair 即 2400 万条记录，1.1 GiB，`CB` / `UB` / `NH`，2 万个基因，`examples/quant_bench.rs`，release 构建，32 核机器，内存预算 16 GiB，不 spill）：

| 产物 | 1 线程 | 4 线程 | 8 线程 | 16 线程 |
|---|---|---|---|---|
| 覆盖度（bigWig） | 15.0 s（扫描 6.1 s） | 5.0 s | 4.9 s | 5.0 s（扫描 2.4 s） |
| 单细胞基因（UMI） | 11.2 s | 3.3 s | 2.7 s | 2.6 s（9.1 M 记录/s） |
| scATAC（片段 + tiles + gene activity） | 28.3 s | 7.6 s | 6.0 s | 5.7 s |
| 以上全部，一次扫描 | 44.3 s | 13.1 s | 10.8 s | 8.9 s |

4 线程以上的瓶颈都在串行部分：覆盖度是合并线程上的事件扫描（约 2 s）和 bigWig 写出（约 2.5 s，`BigWigWriter` 的切块与 zoom 在调用线程上）；基因计数是读取线程切记录（约 10 M 记录/s）；ATAC 是各累加器在合并线程上各自给条形码重新编号。

## 5. 待设计

- Source 接口：catalog、对象依赖、decode kernel 的 Rust trait 与 C ABI
- Planner 的通用对象模型与代价模型
- Output schema：规范行的类型体系与变长输出
- Python 层之上：数据源声明（清单 / 配置）、访问函数构造、与 PyTorch 的衔接
- 派生产物的缓存管理与失效策略

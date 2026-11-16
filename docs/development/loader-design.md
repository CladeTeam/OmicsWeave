# 统一 Loader 设计方向

[文档首页](../README.md) · [当前架构](../design.md) · [训练示例](../../examples/README.md)

**状态：设计提议。** 当前发布的训练入口是普通 PyTorch `DataLoader` 配合格式 reader，见 `examples/`。本文中的 Source catalog、通用 planner、output ring、异步 I/O 编排均不是当前 OmicsWeave 的公共 API。CellPack 模块表和 CENOMultiome / NonB 数据形态是设计时的背景记录，不表示已验证其外部项目当前状态。

## 提议的 N + 1 + M 架构

### 适用场景

CellPack 的 loader 模块虽然随 cpk 格式一起开发，但它解决的问题并不依赖 cpk。它适用于同时满足以下条件的场景：

- **可随机访问**：每个数据单元（item）都能定位到它依赖的字节范围或其他对象；
- **访问序列已知**：要读哪些 item、按什么顺序读，在 run 开始前就确定，或者可以由确定性函数按区间计算（例如按 seed 生成的 shuffle）；
- **可缓存**：重复访问或相邻 item 会共享底层数据，缓存能省下 I/O 或解码开销。

训练数据加载正好满足这三点。由于访问序列已知，I/O 和缓存可以离线或在线（带 lookahead）提前编排，做到 mmap 加 OS 启发式无法做到的预取和淘汰。

如果每种 a 端格式都配一个面向每种 b 端训练格式的加载器，就是 N × M 的组合。OmicsWeave 把加载拆成三部分，得到 N + 1 + M：

| 部分 | 数量 | 职责 | 不负责 |
|---|---|---|---|
| **Source 适配器** | N（每种 a 端格式一个） | catalog：把文件组织进全局 item id 空间，给出 item 依赖的对象（文件字节范围、解码块等）及其尺寸；decode kernel：把读入的字节解码成该 source 的规范行 | 调度、缓存、I/O、输出 |
| **通用运行时** | 1 | Accesses（position → item 的纯函数）；Planner（根据已知序列、对象尺寸和代价模型决定缓存驻留）；Executor（worker 池，显式预取、等待 frame、解码、映射）；Output ring（固定 slot、多 schema region、lease 背压、多 rank 分发） | 任何格式细节、训练语义 |
| **Output schema / process** | M（每种 b 端训练格式一个） | 把规范行转换成固定大小的输出记录：dense 向量、one-hot 序列、token 序列等 | 数据从哪里来 |

### CellPack loader 的现状

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

### 推广时要解决的问题

1. **多级派生对象。** CellPack 有两级缓存对象：原始页和解码行。组学格式需要更多层级，例如 BGZF 的原始页 → 解压块 → 解析出的记录 → item 输入。Planner 的对象需要泛化成「带 key、尺寸和生成代价的对象」组成的依赖链。
2. **item 依赖多个对象，甚至跨多个 source。** FASTA 窗口跨越多个页，相邻窗口互相重叠；一个 consensus item 同时依赖 ref 区间、VCF 区间和常驻的编译结构。
3. **常驻对象。** 有些 source 的索引或编译结构整体常驻内存（fai、VCF 编译结构、基因注释），不参与 frame 调度，但要计入内存预算。
4. **变长输出。** 输出 ring 的每个 schema 使用固定的 `cell_bytes`。变长序列、CSR 行等需要 padding 加长度，或者设计 ragged region。
5. **多组学对齐（Weave）就是 Accesses。** 一个 position 的 S 个 access 可以分别指向不同 source 的 item，CellPack 的多 schema 已经支持每个 schema 取不同的 cell。因此 sample / cell / gene / region 之间的 join 归结为构造访问函数，不需要单独的执行路径。
6. **不可随机访问的格式。** 普通 gzip 无法给出 item 的字节范围，需要派生索引（gzip checkpoint 索引）或整体常驻的解压结果，遵循先转换、再缓存的处理顺序。

CellPack 与 pyconsensus 是已经成熟的独立项目，OmicsWeave 只借鉴它们的设计，不合并、不修改它们的代码。通用运行时在 OmicsWeave 内独立实现。

## 起点用例：CENOMultiome 与 NonB

### 两个项目做什么

- **CENOMultiome**：冻结的 DNA encoder 加多组学 any-to-any 模型，DNA + 已观测组学 → 缺失组学。54 个细胞，hg38，131K（以及 512K / 1M）窗口。
- **NonB**：Non-B DNA（G4、R-loop）。分两条管线：peaks 是 5 kb 窗口二分类；tracks 是 131 kb 窗口的 bigWig 信号回归，最多 146 路 track。数据与训练编排由项目自己的数据层完成。

### a 端：涉及的输入格式

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

### b 端：训练样本的形态

两个项目的训练样本都以**基因组窗口**为单位，同一窗口上叠加多层坐标对齐的数据：

| 层 | CENO | NonB |
|---|---|---|
| 序列 | uint8 ASCII tokens `[L]` + valid，L = 131072（DNA 上下文可扩到 512K / 1M，监督中央 131K） | peaks：one-hot `(L, 4)`，L = 5000；tracks：字符串，L = 131000，送 CENO tokenizer |
| 1 bp 轨道 | 每个细胞最多 9 路（ATAC、DNase、CAGE±、RNA polyA±、RNA total±），各带 valid mask | `(L, n_tracks)` 的 bp 级信号，n_tracks 最多 146 |
| 粗分辨率轨道 | 128 bp 的 histone（11 种 mark）和 CTCF，各带 valid | 训练端 `avg_pool1d` 成 `resolution` bins |
| 二维 | Hi-C 接触矩阵（2048 bp bins）+ mask；1 kb / 5 kb 输入 | — |
| 标签 / 掩码 | observed / target 掩码（任务格点）、dropout 区间、每 track 缩放均值、留一 atlas 参考 | 窗口标签（任一 track 覆盖即为 1）、`(n_tracks, nbins)` 0/1 覆盖 |

### 共同的样本构造逻辑

- **窗口采样**：CENO 用固定锚点区间加 shift 增强；NonB tracks 整染色体平铺、相邻窗口随机重叠、每 epoch 用 `seed + epoch` 重新覆盖；NonB peaks 正负例采样，按 locus 分折。划分按整条染色体或 locus 分组，防止泄漏。
- **逐层变换**：归一化（bigWig 头部 `sumData / nBasesCovered`，或全基因组非零均值）、重复聚合（max / mean / merge / intersect）、背景扣除（`log2((a+1)/(b+1))`）、分 bin 求均值、bf16 量化、shift 时 128 bp bin 的线性插值。
- **反向互补**：序列取互补并翻转，同时交换正负链 track。
- **多细胞共享序列**：CENO 中同一位点的多个细胞共用同一段序列（`sequence_key`）。

### 对 OmicsWeave 的启示

1. **item 是「坐标系上的区间查询」，不是行号。** 访问函数给出 `(contig, start, end, shift, strand, cell/sample)`；不同于 CellPack 每个 cell 对应独立记录，这里相邻窗口、平铺重叠和 shift 会大量共享底层数据，按基因组分块缓存的收益很高。
2. **主流格式自带两级结构，适合「原始块 → 解码块」的缓存模型。** bigWig 是 R-tree 索引加 zlib 压缩数据块，`.hic` 是块索引加压缩块，BGZF / tabix 同理，FASTA 由 fai 直接给出字节偏移。这些格式可以按原样做显式 I/O 和预取，不必先转成自有格式。
3. **CENO 离线转 zarr、NonB 在线读 pyBigWig，是同一个问题的两种妥协。** 前者产生海量分块文件，还要整套清单来管理；后者每个 worker 进程逐窗口读远端存储上的 bigWig，按 float64 逐碱基估算内存，需要大量 CPU 和内存。OmicsWeave 要做的是直接从 bigWig 按已知访问序列预取和缓存，让在线读取够快，从而不需要离线转换。
4. **变换应当是可组合的层操作。** 归一化、重复聚合、背景扣除、分 bin、量化、RC 换链在两个项目里反复出现，属于 M 端（或 a、b 之间）的通用 process，不应由每个项目各写一份。
5. **Hi-C 是二维 item**：一个窗口对应一个 `n × n` 的子矩阵，输出形状和缓存粒度都与一维轨道不同，需要单独设计。

## 尚需设计的接口

- Source trait / C ABI：catalog、对象依赖、decode kernel。
- Planner 的通用对象模型、代价模型与驻留策略。
- 固定及变长输出 schema、lease 和背压。
- 数据清单、访问函数和多进程共享方式。
- 统一派生产物身份和失效机制。

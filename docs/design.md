# 架构与设计边界

[文档首页](README.md) · [Rust 接口](rust-api.md) · [统一 Loader 提议](development/loader-design.md)

## 从文件到训练数组

例如，训练需要同一段 131072 bp 的 DNA 与 ATAC 信号：FASTA 的索引定位序列块，bigWig 的索引定位信号块；reader 解码后返回 `uint8 [131072]` 碱基和 `float32 [131072]` 信号。样本采样、归一化和组装由调用方完成，无需预先展开整个训练集。

当前实现包含格式 reader / writer、格式转换、基因注释和 BAM 定量分析。训练通过 Python reader 与普通 PyTorch `DataLoader` 组合；带通用 Source、planner 和 output ring 的统一运行时仍是[设计提议](development/loader-design.md)。

## 代码分层

| 路径 | 职责 |
|---|---|
| `crates/core` | 坐标、访问 trait、输出布局、错误、I/O 抽象、执行器和缓存 |
| `crates/bgzf` | BGZF、虚拟偏移、BAI / CSI 共用结构 |
| `crates/fasta`、`bbi`、`hic`、`vcf`、`sam`、`bam`、`cram`、`cool`、`bed` | 格式解析、编码及各自支持的查询接口 |
| `crates/htscodecs-sys` | 静态链接 CRAM 所需的 htscodecs |
| `crates/consensus` | BCF 变异应用到参考序列 |
| `crates/convert` | 有线程和内存预算的格式转换、外排序和输出发布 |
| `crates/gtf` | 注释列式表、区间索引、派生特征 |
| `crates/quant` | BAM 一次扫描、多种累加器、矩阵和轨道产物 |
| `crates/python` | PyO3 绑定、NumPy 存储和异常映射 |
| `src/omicsweave` | Python 类型、路径与格式对象、查询、CLI 和任务组装 |
| `refs` | 差分测试的参照实现，`publish = false` |

目录和 crate 包名分别使用短名与 `omicsweave-<name>`。`cool`、`bed` 等组件可服务转换，并不代表 Python 提供同名格式的区域 reader。具体能力见[格式表](formats.md)。

## 设计原则

### 沿用标准格式和索引

原始数据是事实来源。优先使用 FASTA / 2bit、bigWig、bigBed、hic、BAM / CRAM、BCF 和标准索引。普通 gzip FASTA、VCF、SAM、BED、cooler 等先转换成适合目标查询的格式。

`.fai`、`.gzi`、`.bai`、`.csi`、`.crai` 和 `.sidecar` 都是可重建的伴随产物。索引的查找规则见[格式与索引](formats.md)；sidecar 的抽样指纹与查询缓存的进程内 namespace 是两种不同机制，均不应解释为任意源文件改写的自动检测。

### 查询、转换与分析各负其责

- 查询按坐标读取文件已存储的数据；默认不做过滤、归一化或计数。
- 转换改变存储格式，必要时排序；不把 BAM 覆盖度分析混入转换路由。
- 定量层显式选择过滤、特征分配、计数和去重规则，详见[定量分析](quantification.md)。

Consensus 是明确的派生查询：读取参考和 BCF，按选择的样本、单倍型及选项生成序列。参考和变异只加载查询所需部分。

### 显式 I/O 和可组合的内存接口

Rust 格式层通过 `ReadAt` 等接口读取所需范围，保留描述、读取、解码和写入调用方缓冲区的能力。它不以 mmap 作为格式查询接口。当前文件读取器并不因此具备通用异步 I/O planner、`O_DIRECT` 或跨进程 output ring。

固定大小的数值输出可以直接写入调用方存储。变长输出通常先在 reader 工作区暂存，得到大小后写入目标。reader 仍可持有元信息、索引、解压器和可复用工作区；“外部分配输出”不表示全部内存都由调用方分配。

### 缓存所有权明确

固定大小的查询默认分配 NumPy 结果，也可写入调用方的 `out=`。缓存、解码工作区和结果数组拥有各自的生命周期。`clone()` 共享源会话和缓存，使用独立执行状态；独立打开同一路径不会共享 namespace。

缓存是进程内对象，不支持 pickle 或 fork 后复用。DataLoader 每个 worker 创建自己的缓存和 reader，见[缓存指南](query-cache.md)与[训练示例](../examples/README.md)。

### 在线查询与文件产物并存

训练可直接读取原始格式构造样本；转换和定量则显式生成文件。转换与定量使用三段式流水线：串行读取 / 分批、并行计算、按序合并与写出。输出顺序与线程调度分离，外排序用于控制排序缓冲的大小。

`memory` 是引擎对主要常驻数据、在途批次和排序缓冲的预算，不是操作系统的 RSS 限制。缓存的 `max_bytes` 只核算其管理的载荷，也不是整个 reader 的内存上限。

## 兼容性证据

测试调用进程内的参照实现：htslib、samtools、bcftools、libBigWig、lib2bit、straw 和 featureCounts。它们不是用户运行 `ow` 时需要安装的 CLI。运行时所需 native 依赖由 Rust 构建配置静态链接；分发仍需核查每项依赖的许可证。

“兼容”按具体产物定义：有的比较逻辑记录或数组，有的比较索引内容，有的比较未压缩流或完整文件。不能将一种比较的通过推广为所有压缩文件都逐字节相同。已知差异分别记录在[Rust 格式语义](rust-api.md)、[转换](conversion.md)、[注释](annotation.md)和[定量](quantification.md)中。

源码测试是兼容性检查的入口；历史合成基准不是当前 release 的验收结果。验证流程和基准边界见[开发与发布](development.md)及[基准说明](benchmarking.md)。

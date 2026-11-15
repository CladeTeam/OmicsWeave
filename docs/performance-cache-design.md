# 性能优化与外部缓存接口设计草案

状态：方案已确认，第九轮正在实施。本文依据源码审查和 R8 不复用跨查询解压结果的基准；具体落地接口及每包进展见 [执行记录](optimization-round9/coordination.md)。

## 1. 目标及验收边界

优化顺序为：减少算法工作量和重复工作，减少 I/O/解压/解析，减少分配/复制/跨语言调用，最后优化有证据支持的热循环。缓存不能替代无缓存路径的优化。

缓存接口服务于未来统一 loader，只覆盖当前已经提供 query 接口的数据格式。COOL/mcool 等仅参与 convert 的格式不设计缓存接口；其工作仅为转换算法、I/O、并行执行及内存优化。底层存在可重复解压的数据块，本身不构成新增缓存接口的理由。

“所有场景优于对应 binding”作为逐项验收目标：固定输入、输出语义、调用方式、缓存条件和运行环境，任何慢项都单独分析；不以整体平均数掩盖慢项，不用批量接口的结果替代单次调用结果。不能预先保证任意数据、硬件和调用方式都绝对更快，也不能以功能缩减达成表面优势。

现有数据主要是合成数据，仍需增加真实文件、退化输入和不同大小的数据集。性能状态分为领先、未分出稳定差异、落后、无对应接口；后两种不能标记为性能目标已达成。

## 2. 已核实的起点

| 项目 | 当前情况 | 对下一轮的含义 |
|---|---|---|
| BGZF FASTA | reader 复用 raw/text 容量，每次重新解压 | 加入外部解压块复用入口；保留现有连续读取快路 |
| HIC | 保留 zoom/site 元数据，每次 inflate 并遍历块 | 复用解码记录及其查询索引，避免命中后仍扫描整个大块 |
| BBI | 已公开块准备接口，普通 reader 每次读取/解压数据块 | 缓存不可变块和准备元数据，复用已有视图设计 |
| COOL | 已按 HDF5 chunk 显式读取，调用线程持有 scratch | 仅优化转换任务划分、filter/类型转换和临时缓冲，保留 libhdf5 fallback |
| 2bit | N/mask 已二分查询；逐次读取 packed bytes 并解码 | 先分解窄窗固定开销，不把已有二分再列为新增优化 |
| Python | 多个查询在 Python 分配 NumPy、调用 native、构造结果对象；native 状态经 Mutex 并释放 GIL | 微秒级查询应实测这些固定成本，缓存本身不能解决它们 |

最近基准中的优先慢项：FASTA 128bp 为 5.223µs，对 pyfaidx 3.229µs；2bit 128bp 为 5.317µs，对 py2bit 0.904µs；2bit 4096bp 为 6.253µs，对 py2bit 4.409µs。这些是定位入口，不是对瓶颈原因的测量结论。

出处：[查询报告](optimization-round8/query-performance.md)、[FASTA reader](../crates/fasta/src/reader.rs)、[2bit 解码](../crates/fasta/src/twobit.rs)、[HIC reader](../crates/hic/src/reader.rs)、[BBI reader](../crates/bbi/src/reader.rs)、[Python 共享状态](../crates/python/src/shared.rs)。

## 3. 所有权和职责

沿用 [总体设计](design.md) 的显式 I/O、描述/操作分离及外部分配原则：

1. 格式组件计算块依赖、说明大小/布局、验证并解码到调用方提供的 frame，从 frame 视图向调用方输出写入数据。
2. 外部 cache/planner 分配和持有 frame，决定准入、预算、预取、共享及淘汰。格式组件不创建私有 LRU，不持有全局缓存，不在后台自行预取。
3. reader 的便捷接口可以连接调用方显式传入的缓存对象；reader 只保留访问句柄，不取得缓存策略的控制权。关闭一个 reader 不清空外部共享缓存。
4. 已有 reader 的临时工作缓冲与解压器继续作为执行状态；它们的容量复用与数据缓存分别计量。底层描述/解码接口仍可由调用方直接组合。
5. 缓存中存放独立于查询窗口的规范数据。首次实现不缓存完整 Sequence、Signal、ContactMatrix 或按窗口裁剪的最终结果。

命中路径：描述依赖 → 获取已验证块的只读 lease → 范围定位 → 写入输出。缺失路径：外部分配/调度读取 → 解压/解析/准备 → 验证成功后发布 → 写入输出。

同一次查询内，一个块只解压/准备一次；相邻区域或不同 reader 可借助同一外部缓存复用这个块。

## 4. 哪些接口支持缓存

| 包/输入 | 缓存单元与内容 | 决策 |
|---|---|---|
| fasta / BGZF | 物理 BGZF member 的解压文本及必要的映射元数据 | 支持；不能把一个窗口误当成块 |
| fasta / plain FA | 已有 FAI + 显式范围读取，处理换行后写输出 | 不增加专用数据缓存接口 |
| fasta / 2bit | 可选固定范围的已解码 tile；大写规范碱基，N 语义保留，soft mask 在输出时应用 | 支持可选 tile 复用；tile 大小/准入由外部决定，无缓存窄窗优化另做 |
| hic | 版本/染色体对/unit/resolution 下的逻辑块，已解码记录及查询索引 | 支持；同一内容可服务 matrix/pixels，转置和裁剪在输出阶段做 |
| bbi / bigWig | 块数据与已验证的 WigBlock/ZoomBlock 查询元数据 | 支持；缓存命中也不应反复全块 prepare |
| bbi / bigBed | 块文本、记录偏移及区间查询元数据 | 支持；未压缩 bigBed 也可能因重复解析受益 |
| bam / indexed BAM | BGZF member 解压文本，BAI/CSI 选择的 chunk 共享底层块 | 新增 AlignmentReader 已有 query，支持外部块缓存及独立上下文 |
| vcf / indexed BCF | BGZF member 解压文本，CSI 区域及样本选择后的 GT 输出 | 新增 VariantReader 已有 query，复用 indexed BGZF 缓存，不缓存最终矩阵 |
| cool / mcool | 仅为 convert 输入 | 不新增缓存接口；优化 chunk 转换与执行划分 |
| bgzf | FASTA 等查询路径使用的通用 codec | 优化 codec 和外部分配操作，不给仅用于转换的顺序流新增缓存接口 |
| vcf / 顺序流 | 顺序 VariantReader 用于转换及装载 | 顺序路径不强加缓存；indexed query 单独按上面的协议接入 |
| sam | 当前职责为 SAM → BAM 编码 | 不增加没有重读场景的数据缓存；复用编码器/scratch |
| consensus | 按区间经 BcfReader 读取记录；引用的序列通过 SequenceRead 获取 | 继承 reference 的缓存能力；batch 同区域共享一次记录解码，不缓存最终各样本序列 |
| convert | 默认顺序流转换 | 不默认污染随机访问缓存；共享 codec/scratch、检查额外全表扫描 |
| core / python | 公共描述/lease 协议与外部控制入口 | core 不依赖格式包；Python 不在每个块上回调 Python dict |

先确认格式已有 query 接口，再判断重复 I/O、解压、解码或准备是否值得缓存。普通 FA 的范围读取不需要强加一层小对象缓存，convert-only 格式不纳入这套接口。

## 5. 公共协议

### 5.1 最小概念

- `SourceToken`：一次已打开、按不可变内容使用的数据源会话的身份。与文件路径、抽样 hash 区分。
- `BlockKey`：source token、格式/解码布局版本、格式内部流或 dataset 标识、块身份及必要的解码配置。使用结构化字段，不拼接字符串。
- `BlockRequest`：key、需要读取的物理范围、已知的解码长度或未知标记、对齐、准备元数据需求及资源限制。
- `BlockLease`：外部持有 frame 的只读借用/引用计数句柄。lease 活跃时内容和准备元数据不可修改，淘汰不能使其失效。
- `DecodeScratch`：每执行上下文独立的临时缓冲和解压器，不算作可复用数据块。

这些名称是讨论用名称，不代表承诺新增同名公开类型。先做 BGZF 和 HIC 两个原型检验生命周期，再冻结 Rust 签名。

### 5.2 格式层需要能表达的操作

```text
describe(query, caller_request_buffer) -> 块依赖及输出布局
decode(request, stored_bytes, caller_frame, caller_scratch) -> 已验证的准备元数据
write(query, leased_block_views, caller_output) -> 完整结果
```

这是语义示意，不是需要所有包照抄的 trait。优先补齐和组合查询路径现有 Gzi::blocks、ZoomData::blocks_for、BBI prepare 能力，避免另建一套重复解析器。

请求列表也允许复用调用方缓冲。plan 只持有稳定 key、范围和数值元数据，不持有易因 reader 下一次调用而失效的裸指针。未启用缓存时继续调用直接读取/解码路径，避免查表、Arc、互斥及对象分配。

### 5.3 缓存层需要的操作

| 操作 | 语义 |
|---|---|
| lookup / acquire | 命中返回只读 lease；缺失不伪造空块 |
| reserve / publish | 外部先做预算准入和 frame 分配；完整解码/验证成功后原子发布 |
| evict(key) | 从可命中集合移除；活跃 lease 仍有效 |
| clear(namespace) | 清除指定源或全部条目，同时阻止早于 clear 的在途加载重新发布 |
| resize(bytes, entries) | 显式调整预算；说明 pinned 条目如何处理，不隐式超额 |
| stats | 查看命中、缺失、准入拒绝、淘汰、实际读取/解压次数、驻留/pinned/在途内存 |

任意 safe Rust 缓存实现不能通过返回错误 key、过短 frame 或别名可变内存触发格式组件 UB。准备结果必须和其不可变 bytes 绑定；构造入口受检，元数据用偏移而非自引用指针，取回时验证 key/layout 一致。缓存命中不等于可以信任任意外部字节。

可提供一个显式创建的内存缓存便于 Python 使用；Rust 端允许其他外部实现。默认无缓存，首版不引入磁盘缓存、跨进程共享内存或持久化已准备块。

## 6. Python 建议接口

以下仅为调用外形；新增对象由调用方创建，容量/准入/淘汰策略均属于这个对象。

```python
cache = ow.cache.BlockCache(max_bytes=256 << 20, max_entries=8192)
reader = ow.query.SequenceReader("ref.fa.gz", cache=cache)
contacts = ow.query.ContactReader("map.hic", cache=cache)

sequence = reader.fetch(region)
matrix = contacts.matrix(window)

cache.stats()
cache.evict(key)                  # key 来自计划/条目描述接口
cache.clear(source=reader.cache_namespace)
cache.resize(max_bytes=128 << 20, max_entries=4096)

reader.close()                    # 不清空其他 reader 使用的 cache
cache.clear()
```

- 同一个外部 cache 可以对多个格式和文件统一预算，但各格式条目按 namespace/layout 隔离。
- plain FA 等没有对应缓存能力的输入，显式传入数据缓存时报清楚的参数错误，不静默假装生效。
- 底层计划/块接口允许 loader 自己选择预取和发布时机。首版可先开放同步组合，不启动隐式预取线程；普通 fetch 在缺失时同步加载。
- 热路径完全在 Rust 执行；Python 控制方法用于配置和管理，不逐块调用 Python 回调。
- 为缓存共享增加 reader 执行上下文克隆能力：共享文件会话/不可变索引和 cache，scratch 独立。多个线程共享当前同一个 reader 仍可串行使用，不能仅凭“释放 GIL”宣称并行。
- DataLoader pickle/spawn 仅保存重新打开的配置，不序列化缓存内容或进程内 source token；子进程显式创建/绑定自己的 cache。OS fork 后不复用可能继承锁状态的 reader/cache。
- 保持现有返回类型、排序/裁剪语义及 out 用法；不能让普通返回值因 cache 被淘汰或下一次查询而改变。

## 7. 身份、失效与正确性

缓存会话要求输入文件和对应索引在会话期间保持不变。任意 ReadAt 本身没有不可变保证，不能默认给所有 ReadAt 开启长期缓存。

默认每次独立打开生成新 source token；同一个源会话的显式克隆可以共享 token。不同文件、不同版本、不同 reader 的独立打开不会仅凭路径/长度相同而复用。需要数据更新时关闭并重新打开 reader，重新加载元信息和索引；仅清空块缓存不能修复陈旧索引。

当前 sidecar SourceIdentity 是抽样指纹，源码已说明等长局部改写可能漏检，不能将它升级为严格内容身份。mtime/inode 也不是不可变内容保证。不在每次热查询中重算文件 hash；第一阶段不支持对正在原地改写的文件进行一致快照读取。

格式细节必须保留：

- BGZF：CRC/ISIZE 与文本映射验证；GZI 可略过空 member，物理块身份与逻辑文本位置要区分。缓存的是完整解压文本，处理 FASTA 换行仍由 decoder 完成。
- HIC：块身份包括染色体对、unit、分辨率及版本；缓存源坐标记录，不缓存查询方向。物理反序、重复物理引用、重复像素和原有矩阵赋值规则保持。索引/排序不得改变重复记录的既有语义。
- BBI：保留输入顺序、嵌套区间、重复起点及 last-write-wins；未压缩块也必须验证边界。准备元数据与 bytes 生命周期绑定。
- 2bit：tile 粒度固定在同一缓存命名空间中；不同 tile 布局不能共享 key。N/soft mask 行为与直接解码一致。

## 8. 内存和并发

1. 同时限制字节预算和条目数，空块、元数据和大量小块也受限。计费使用分配容量及准备索引占用，不只用有效 payload 长度。
2. 分别统计可驻留字节、活跃 lease 对应字节、在途解码预算及 scratch；不要将 cache 的 max_bytes 宣称为进程峰值 RSS。
3. 缓存支持的 frame 被淘汰后，只要有 lease 仍需计费。预算不足、超大单块或无可淘汰项时拒绝准入，正常走有界 scratch 路径；不能挂住查询等待自己持有的 lease。
4. 解码尺寸未知的 HIC 等请求不能伪造精确大小。使用外部解码预算、可增长受限 frame 或 sidecar 的尺寸提示；发布前核对实际长度。sidecar 不替代内容验证。
5. 查表/发布的锁不覆盖读取、解压和输出拷贝。相同 key 的并发缺失合并一次加载；失败会唤醒等待者，不缓存半成品，不让后续查询永久失效。
6. clear 使用 generation 隔离在途提交；旧 lease 可完成读，新查询不能命中过期条目。需要显式说明 clear 后暂时保留的 pinned 内存。
7. 首版只提供最小可测的外部缓存实现；具体 LRU/CLOCK/已知访问序列策略可替换。预先知道访问序列的 loader 应能直接控制准入/淘汰，不被固定 LRU 限制。

## 9. 复杂度级别的优化目标

| 方向 | 优先检查的重复工作 | 目标与限制 |
|---|---|---|
| BGZF | 跨查询重复读取/解压同一块 | 命中后只做范围定位和 O(L) 输出；无缓存和 miss 保留合并读取/直接输出，避免一块一次系统调用 |
| HIC | 命中后全块解码/扫描，pixels 每次排序 | 一次建立记录索引，查询访问相关行/候选记录；有序块结果可研究合并代替每次全量排序。分别报告首次准备、每查询候选数和输出量，不无条件承诺 O(log N + K) |
| BBI | 重复 prepare，大块内全扫描，重叠区间多次覆盖输出 | 复用已准备范围索引；按排序与嵌套实际结构选路径。区分真实命中 K 和保守候选 C，保留已存在的复杂度改进 |
| COOL 转换 | 任务边界导致重复 filter/decode 同一 chunk、重复数值转换 | 按 chunk 组织转换任务、保留三列独立边界，复用线程 scratch；不新增缓存接口 |
| Consensus | 多 selector 反复取参考、重复选位点/解析 selector | 同区域共享加载和准备；输出长度仍有 O(L) 下界，不缓存最终完整样本序列掩盖算法开销 |
| 转换/编码 | 全文件扫描、多轮排序、重复字符串解析和中间全量副本 | 优先消除额外整轮工作，保持流式、有界内存及既有输出语义 |

L 为输出规模。输出矩阵的初始化/写入不可简单省略，除非 API 明确接受稀疏输出；稀疏 API 的性能不能代替密集 API 的验收。

## 10. 系统性定位窄窗性能

把一次查询拆为参数处理、Python/FFI、锁/GIL、输出分配、范围定位、读取、解压、块内查询、输出写入和结果包装。分别测组件成本与端到端成本，组件计时只作解释，不混入对用户的速度比表。

优先原型：默认公开查询在 native 内完成最终 NumPy 存储分配及填充，减少分开的 Python 分配/native 调用；研究紧凑结果包装、减少重复字段提取/名称解析/跨层校验，以及小查询释放 GIL 的阈值。是否采用由实测决定，不能在持有 GIL 时等待一个需要 GIL 才能完成的线程释放共享锁。任何阈值必须验证多线程、close 和异常行为。

2bit 的 lookup 解码和 N/mask 已有优化基础，需证明其在当前窄窗中占比后再改 SIMD/unchecked。新增缓存不会被当成消除这些固定成本的替代品。

## 11. 基准与接受标准

保留两个独立问题：

- 解码/查询主表：延续已加载 reader、无跨查询解压结果命中的口径；每个场景只有一个我们的结果，对应 binding 及 CLI 写文件后读取。
- 缓存专项：外部显式启用缓存，列首次缺失、稳定命中、相邻窗口复用、容量不足/抖动、单次顺序扫描。对照缓存行为和预算均披露，能控制时对齐；不能对齐时不宣称完全同条件。

每种条件分别列单值，不合并成跨格式、跨窗口或不同绑定的最小–最大范围。缓存专项至少报告实际物理读取字节、解压次数、命中率、准备次数、分配数、预算和峰值内存，证明收益来自哪里。

扩充空/短/长窗口、块边界、稀疏/密集、嵌套/重复、随机/顺序/热点/无复用、不同 cache 容量、单线程/多线程、真实文件和合成对抗输入。正向及 transposed HIC、FRAG 等查询功能不可通过少测而隐藏。COOL fallback 等转换路径另行验收。

每个慢项形成记录：准确输入/版本、端到端差距、主要耗时层、优化候选、修改后的同条件结果和剩余差距。共享 CPU 测量不确定时做交错重复和对照，而不是把噪声范围当作收益。发生崩溃、输出不一致或语义降级直接拒收；未达性能目标则保持未完成状态。

## 12. 后续每包一个 subagent 的执行顺序

按协议依赖分批安排任务：先验证公共缓存协议，再接入各格式并统一验收。

1. 主代理与 core 负责人先完成最小协议；fasta/BGZF 与 hic 做两种原型验证 frame/lease 生命周期及无缓存分支成本，避免所有包同时发明接口。
2. 12 个包分别安排一个新负责人，按并发槽位分批；每人只修改分配的包。共享 core、Python 绑定和基准脚本由对应 owner 统一改，跨包需求通过协议协作。
3. bbi 依赖公共协议冻结后接入；cool 仅优化转换路径，bed/vcf/sam/consensus/convert 开展各自算法和工程优化，不给 convert-only 格式强加缓存。
4. 每包提交复杂度变化、正确性边界、unsafe 不变量、内存/缓存统计及配对性能证据。主代理统一端到端验收，逐个慢项收敛。
5. 编译、测试、基准和重任务仅在届时动态发现并验证的现有 CPU pod 执行，不在开发机执行、不自行申请 pod。委派时传递实际 SSH、realpath 及这些约束；不复用过期 pod 地址。

新增 unsafe 仅限有测量依据的热循环，配套 SAFETY 说明；输入长度、布局、算术、缓存键/版本和元数据保持 safe 验证。不能用缓存命中或外部实现提供了某对象，作为 unchecked 安全性的唯一证明。

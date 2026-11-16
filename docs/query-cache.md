# 查询块缓存

[文档首页](README.md) · [区域查询](queries.md) · [缓存扩展方向](development/cache-design.md)

缓存由调用方显式创建和管理，用于让重叠窗口、重复访问及多个查询执行上下文共享底层块。reader 默认不启用缓存，也不自动启动预取线程。

已接入查询格式：BGZF FASTA、2bit、HIC、bigWig、bigBed、BAM、CRAM 与 BCF indexed 查询。普通 FA 直接范围读取。BED、COOL/mcool 及其他 convert-only 格式不提供查询路径，也不提供这套缓存接口。

GTF/GFF3 的 `Annotation` 一次完整加载为只读列式表，外部通过 `Arc<Annotation>` 共享表、索引和派生视图，各自创建 sweep cursor；查询不会重新解压源文件，因此不附加块缓存接口。Rust `IntervalSet` 允许直接编辑公开列，编辑后需调用 `clear_index()`，下次查询再构建索引；它不会自动检测任意列修改。quant 保留自身参考窗口与任务内存预算，不套用格式 reader 的缓存控制面。

## Python：挂上缓存、克隆与计数

Python 只负责顶层用法：创建 `BlockCache`、把它交给 reader（`cache=` 或 `set_cache`）、用 `clone()` 给其他线程独立的执行状态、读 `stats()`，以及 `clear()` / `resize()`。哪些块被某次查询依赖、何时预取和淘汰，是 Rust 层（未来 loader）的接口：`describe_cache`、`prefetch_cache`、`BlockKey`、按源会话清理、单条淘汰和读取计数都只在各格式 crate 中提供，不暴露到 Python。

```python
import omicsweave as ow

cache = ow.cache.BlockCache(max_bytes=256 << 20, max_entries=8192)

with ow.query.SequenceReader("reference.fa.gz", cache=cache) as reader:
    sequence = reader.fetch(ow.Region("chr1", 1000, 2000))  # 缺失的块同步加载并准入
    clone = reader.clone()

# 关闭原 reader 后，clone 仍保有同一源会话，使用独立的工作缓冲。
try:
    next_sequence = clone.fetch(ow.Region("chr1", 1500, 2500))  # 命中重叠的块
    print(cache.stats())
finally:
    clone.close()

cache.clear()
```

命中复用的是完整解压块或准备好的记录 / 索引，不是完整窗口结果。`TrackReader`、`StrandedTrackReader`、bigBed `IntervalReader`、`ContactReader`、`AlignmentReader`（BAM 与 CRAM）、`VariantReader` 和 `ConsensusReader` 都接受 `cache=`、`set_cache`；同一 cache 可同时服务多个 reader。CRAM 的缓存保留不可变的容器元数据、压缩数据及解码 slice，不同字段选择读到的记录一致；reference 与文件、索引一起属于不可变源会话，clone 共享会话，换文件或 reference 则重新打开。

ConsensusReader 使用已打开的 SequenceReader 作为 reference 即可共享它的块缓存；Consensus 按区间经 BCF + CSI 读取记录，克隆共享 BCF 的 header、索引与选项，并使用独立的 BCF 与 reference 执行上下文，不缓存各个样本的最终 consensus。

每个公开的 consensus `fetch` / `fetch_batch` 调用重新加载区域；只有同一次 `fetch_batch` 的多个 allele selector 共享一次 reference 和记录解码。不会因重复请求同一区域而私下命中上次的完整解压窗口。

## 驻留与预算

```python
cache.clear()                                       # 清除全部条目
cache.resize(max_bytes=128 << 20, max_entries=4096)
# 在 reader 打开期间可用 reader.set_cache(None) 显式脱离缓存。
```

便利实现使用 FIFO，命中不更新淘汰顺序。按未来访问序列控制驻留（单条淘汰、按源会话清理）属于 Rust 层 planner 的接口；格式包不会另外维护私有 LRU。

预算同时限制字节和条目数。字节计费包含已准入 payload 容量、准备索引和 entry 值，不包括解码 scratch、准入前在途分配、map/allocator 开销，因此不是进程 RSS 上限。超大条目或预算不足时，查询仍可通过未准入数据或普通路径完成。

`stats()` 区分：

- `resident_bytes` / `resident_entries`：当前可命中的条目。
- `live_bytes` / `live_entries`：包括被淘汰但仍由查询 lease 持有的数据。
- `in_flight`：仍在执行的加载，包括被 clear 隔离的旧加载。
- `hits`、`misses`、`loads`、`waits`、`admission_rejections`、`evictions`：查表与加载控制计数。

clear 不使正在查询的只读 lease 失效；旧加载完成后不会重新发布到已清理的缓存。关闭一个 reader 不清空共享 cache。普通返回结果独立于缓存，后续清理不会修改已返回的数组。

## 源内容与多进程

源文件和索引在会话期间必须保持不变。更新文件后重新打开 reader，重新加载索引；仅清空数据块不能更新陈旧索引。缓存按源会话区分，不按内容 hash，不自动检测任意原地改写。

缓存和打开的 reader 是进程内对象。DataLoader spawn/pickle 只传递 reader 的重新打开配置，子进程显式创建自己的 cache；缓存内容和源会话身份不跨进程序列化。正常 CPython fork 后，继承的 reader/cache 在碰锁前拒绝查询，需在子进程重新创建。直接使用外部 libc fork 绕过 CPython fork hooks 不在此保证内。

POSIX 的原生 child hook 会在用户 Python child callback 之前更新进程标识，包括用户在导入包之前注册的回调。GTF Annotation 与 IntervalSet 的 native 懒索引查询也属于进程内状态；fork 后分别重新读取/反序列化 Annotation，或重新构建 IntervalSet。已经导出的普通列数组仍可作为数组读取。

## 与未来统一 loader 的边界

Rust 格式包提供外部缓存控制、依赖 key、同步预取、独立执行上下文，以及块定位、解码与 prepared view 等低层接口，HIC 另提供具体 block request / 加载接口；Python 只暴露挂缓存、克隆和计数。

本实现尚不构成统一训练 loader，也没有任意异步 frame 向所有 reader 发布的通用协议、可替换 allocator、跨进程共享缓存或磁盘缓存。未来 loader 可在当前控制面和低层操作上接续设计，不应将这些能力视为已经实现。

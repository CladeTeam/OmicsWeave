# 查询块缓存

缓存由调用方显式创建和管理，用于让重叠窗口、重复访问及多个查询执行上下文共享底层块。reader 默认不启用缓存，也不自动启动预取线程。

已接入查询格式：BGZF FASTA、2bit、HIC、bigWig、bigBed、BAM 与 BCF indexed 查询。普通 FA 直接范围读取。BED、COOL/mcool 及其他 convert-only 格式不提供查询路径，也不提供这套缓存接口。

GTF/GFF3 的 `Annotation` 一次完整加载为只读列式表，外部通过 `Arc<Annotation>` 共享表、索引和派生视图，各自创建 sweep cursor；查询不会重新解压源文件，因此不附加块缓存接口。Rust `IntervalSet` 允许直接编辑公开列，编辑后需调用 `clear_index()`，下次查询再构建索引；它不会自动检测任意列修改。quant 保留自身参考窗口与任务内存预算，不套用格式 reader 的缓存控制面。

## 显式创建、预取与查询

```python
import omicsweave as ow

cache = ow.cache.BlockCache(max_bytes=256 << 20, max_entries=8192)
region = ow.Region("chr1", 1000, 2000)

with ow.query.SequenceReader("reference.fa.gz", cache=cache) as reader:
    keys = reader.describe_cache(region)
    reader.prefetch_cache(region)
    sequence = reader.fetch(region)
    clone = reader.clone()

# 关闭原 reader 后，clone 仍保有同一源会话，使用独立的工作缓冲。
try:
    next_sequence = clone.fetch(ow.Region("chr1", 1500, 2500))
    print(cache.stats())
    cache.clear(source=clone.cache_namespace)
finally:
    clone.close()

cache.clear()
```

预取只加载依赖块，不构造最终 Sequence/Signal/Matrix 输出。普通查询在缺失时同步加载。命中复用的是完整解压块或准备好的记录/索引，不是完整窗口结果。

`TrackReader`、bigBed `IntervalReader`、`ContactReader`、`AlignmentReader` 和 `VariantReader` 使用同样的 `cache=`、`set_cache`、`cache_namespace`、`describe_cache`、`prefetch_cache` 和 `clone` 控制面；Contact 的描述/预取参数为其查询使用的 ContactWindow。BAM/BCF 还可通过 `stream_stats()`、`reset_stream_stats()`核对真实读取/解压计数。

ConsensusReader 使用已打开的 SequenceReader 作为 reference 即可共享它的块缓存；Consensus 按区间经 BCF + CSI 读取记录，克隆共享 BCF 的 header、索引与选项，并使用独立的 BCF 与 reference 执行上下文，不缓存各个样本的最终 consensus。

## 驻留与预算

```python
cache.keys()                              # 当前可命中的 opaque BlockKey
cache.evict(keys[0])                       # 单条淘汰
cache.clear(source=namespace)             # 清除指定源会话
cache.clear()                             # 清除全部条目
cache.resize(max_bytes=128 << 20, max_entries=4096)
reader.set_cache(None)                    # 显式脱离缓存
```

key 可以比较、哈希，具有 `source` 属性；不通过手工构造 key 将其他文件数据混入。独立打开的 reader 获得不同 namespace，同一源会话的 `clone()` 才共享 namespace。

便利实现使用 FIFO，命中不更新淘汰顺序。外部 planner 可以根据未来访问序列，通过 keys/evict/clear/resize 控制驻留；格式包不会另外维护私有 LRU。

预算同时限制字节和条目数。字节计费包含已准入 payload 容量、准备索引和 entry 值，不包括解码 scratch、准入前在途分配、map/allocator 开销，因此不是进程 RSS 上限。超大条目或预算不足时，查询仍可通过未准入数据或普通路径完成。

`stats()` 区分：

- `resident_bytes` / `resident_entries`：当前可命中的条目。
- `live_bytes` / `live_entries`：包括被淘汰但仍由查询 lease 持有的数据。
- `in_flight`：仍在执行的加载，包括被 clear 隔离的旧加载。
- `hits`、`misses`、`loads`、`waits`、`admission_rejections`、`evictions`：查表与加载控制计数。

clear/evict 不使正在查询的只读 lease 失效；旧加载完成后不会重新发布到已清理的缓存。关闭一个 reader 不清空共享 cache。普通返回结果独立于缓存，后续清理不会修改已返回的数组。

## 源内容与多进程

源文件和索引在会话期间必须保持不变。更新文件后重新打开 reader，重新加载索引；仅清空数据块不能更新陈旧索引。namespace 是会话身份，不是内容 hash，缓存不自动检测任意原地改写。

缓存和打开的 reader 是进程内对象。DataLoader spawn/pickle 只传递 reader 的重新打开配置，子进程显式创建自己的 cache；缓存内容和 namespace 不跨进程序列化。正常 CPython fork 后，继承的 reader/cache 在碰锁前拒绝查询，需在子进程重新创建。直接使用外部 libc fork 绕过 CPython fork hooks 不在此保证内。

## 与未来统一 loader 的边界

现有接口提供外部缓存控制、依赖 key、同步预取和独立执行上下文。Rust 格式包保留块定位、解码及 prepared view 等低层接口，HIC 另提供具体 block request/加载接口。

本实现尚不构成统一训练 loader，也没有任意异步 frame 向所有 reader 发布的通用协议、可替换 allocator、跨进程共享缓存或磁盘缓存。未来 loader 可在当前控制面和低层操作上接续设计，不应将这些能力视为已经实现。

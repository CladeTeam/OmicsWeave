# 基准测试与结果解释

[文档首页](README.md) · [脚本使用说明](../scripts/README.md) · [开发与发布](development.md)

仓库提供三种互补测量。输入文件、对照依赖、CLI 和运行命令见 [scripts/README.md](../scripts/README.md)。Fixture 和外部工具需要自行准备，完整 fixture 生成器不随当前脚本发布。

## 测量口径

| 入口 / 模式 | 计时包含 | 计时排除 |
|---|---|---|
| `benchmark_simple_queries.py` | 常驻 reader 的一次公开查询和默认结果分配；CLI 启动、写文件、等退出、读回字节 | 导入、打开 reader / 建立内存索引、参数构造、额外 canonical 转换 |
| `benchmark_queries.py`：`raw_reuse` | 常驻 reader 的原生完整结果，包含适配器明确声明的物化工作 | 打开、canonical 转换 |
| 同上：`canonical_reuse` | 常驻 reader 查询加完整输出标准化 | 打开 |
| 同上：`canonical_open_query_close` | 每次打开、查询、标准化、关闭 | 包导入、fixture 准备 |
| `benchmark_query_metadata.py` | 已打开对象的元信息查询 | 打开和数据载荷解码 |

原生结果类型不同，只有完整逻辑输出经过标准化并相等时，才能解释为完成了相同工作。两个查询 runner 在计时前验证 payload，并在每次采样后验证最终输出；失败记录到 `errors`，进程返回非零。缺失 fixture、依赖或 CLI 等初始化错误可能在写 JSON 前终止。

CLI 是真实工具，其成本包含进程启动和输出物化；不能把这类比较称为纯解码器速度。Consensus 的 CLI / pysam 路径按查询打开输入，而 OmicsWeave 使用常驻 reader；报告应保留这个差异。

## 缓存协议

默认依赖热 OS 页缓存，不清全局页缓存。`simple` 的无显式缓存协议对 pyfaidx / Biopython 关闭或限制连续查询间块复用，并检查相邻 BGZF 查询的压缩块互不重叠。

`--ours-cache-bytes` 是另一种协议：OmicsWeave 显式预热块缓存，对照 binding 保留各自默认缓存。两边策略不同，必须同时报告，不能将结果描述为相同缓存实现的比较。plain FASTA、consensus 等不适用的 case 不参与该缓存模式。

探针只诊断某些动态 zlib 调用；它必须在独立的 `--validate-only` 进程中运行，不能 preload 到正式计时进程。未观测到 inflate 调用不单独证明没有解压。

## 应保存的证据

- 待测源码 commit、dirty diff、native 扩展路径与 hash，以及 runner / adapter 的版本。
- Python 与对照包版本、CLI 版本和路径、编译选项、机器与文件系统信息。
- Fixture 的生成规则、seed、文件尺寸 / hash、所有实际查询和 payload 校验值。
- CPU affinity、并发负载、重复次数、每次 wall / CPU 样本与离散程度。
- 缓存协议、冷启动 / 预热、分配模式、CLI 输出目录及文件系统。

JSON 自动记录其中一部分，其余应随结果补充。CPU 时间包含已等待的子进程；固定 affinity 不代表 CPU 独占。输出目录在共享文件系统与本地盘上时，CLI 写入成本可能不同；脚本不执行 fsync，因此不测量持久化提交延迟。

## 结论边界

比较前先检查 `completed`、`errors`、case 的后端状态和 checksum；FRAG 与酶切位点 case 仅做正确性验证。空匹配、缺失后端或错误不能写成通过。

这些脚本覆盖固定合成查询，不测完整 DataLoader、远端对象存储、训练 GPU 吞吐或所有格式；目前没有 CRAM 查询对照组。元信息方法也不能代表读盘和解码性能。历史优化轮次的数据只对其记录的源码、fixture 和环境有效，本页不承诺当前 release 的速度或加速比。

# OmicsWeave 文档

OmicsWeave 按基因组区域读取标准组学文件并返回 NumPy 数组，同时提供格式转换、基因注释和 BAM 定量工具。当前 Python 包与 Rust workspace 版本均为 `0.1.0`。

## 开始使用

| 文档 | 内容 |
|---|---|
| [安装与快速开始](getting-started.md) | 构建要求、安装、最小序列和信号查询 |
| [格式与索引](formats.md) | 哪些格式可以直接查询、需要哪些伴随文件、writer |
| [区域查询](queries.md) | 坐标、链、数组形状、生命周期、线程和进程 |
| [查询块缓存](query-cache.md) | 创建、预取、克隆、预算和失效边界 |
| [格式转换](conversion.md) | `ow convert`、Python API、路由、内存和兼容差异 |
| [基因注释](annotation.md) | GTF / GFF3、区间索引、转录本及派生区域 |
| [BAM 定量分析](quantification.md) | `ow quant`、产物布局、计数定义及参照工具差异 |
| [训练示例](../examples/README.md) | 普通 PyTorch DataLoader 的完整示例 |

## 开发与验证

| 文档 | 内容 |
|---|---|
| [架构与设计边界](design.md) | 当前实现、组件职责和内存模型 |
| [Rust 格式接口](rust-api.md) | trait、布局、格式存储与兼容语义 |
| [基准测试](benchmarking.md) | 脚本入口、输入契约、计时口径和结果解释 |
| [CI 与 PyPI 发布](ci.md) | 平台矩阵、workflow 触发条件和上传流程 |
| [开发与发布](development.md) | 构建、差分测试、打包和发布前检查 |
| [缓存扩展方向](development/cache-design.md) | 已有能力与尚未实现的扩展协议 |
| [统一 Loader 提议](development/loader-design.md) | N + 1 + M 架构及设计背景，尚未实现 |

本文档中的业务文件名是示例，需替换为自己的数据。命令行选项以当前安装的 `ow convert --help` / `ow quant --help` 为准。基准脚本是仓库开发工具，不属于安装后 `ow` 的命令。

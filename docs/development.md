# 开发与发布

[文档首页](README.md) · [CI 与 PyPI 发布](ci.md) · [基准测试](benchmarking.md)

## 构建环境

项目使用 Python 3.12+、Rust 1.89+（edition 2024）、C 编译器和 CMake 3.26+。Python 环境由 `pyproject.toml` / `uv.lock` 管理，Rust 依赖由 `Cargo.toml` / `Cargo.lock` 管理。初次安装会编译 Rust、C 和静态 HDF5。

以下命令在仓库根目录执行。全量构建和测试应在资源充足的计算节点上运行；PJLab 环境使用用户已有且已验证的 worker 与共享盘绝对路径，不在开发盒执行大型编译，也不自行申请 pod。

```sh
uv sync --locked
# Rust 改动后重新构建扩展到当前开发环境
uv run maturin develop --release --skip-install
```

`--skip-install` 配合已有开发环境中的安装使用。切换环境或构建方式后，检查实际导入的扩展，避免把旧二进制的结果算到新源码上：

```sh
uv run python -c 'import omicsweave as ow; print(ow.__version__, ow._native.__file__)'
```

## 验证

```sh
uv run ruff check src scripts .github/scripts
uv run ruff format --check src scripts .github/scripts
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
RUST_TEST_THREADS=1 cargo test --workspace --locked
uv run pytest
```

C reference oracle 的部分测试涉及进程全局状态，串行 Rust 测试有助于稳定对照。`refs/` 中的上游源码只用于测试；可按变更范围先运行目标 crate，再运行发布所需的完整检查。被标为 `ignore` 的大数据测试或手工基准不等于已经运行。

Python / Rust 的关键契约包括：输入和索引校验、空区间与边界、正负链、定长与变长输出、`out=` 的 shape / dtype / 别名、并发与进程生命周期、缓存计费和清理，以及不同线程数和 spill 条件下的确定性。不要把历史通过数当作当前 commit 的验收记录。

## 构建 wheel 和源码包

```sh
uv run maturin build --release --locked --out dist
uv run maturin sdist --out dist
```

应在仓库之外的独立环境安装并测试 wheel，确保导入不依赖 checkout 的 `src/` 或 editable 安装。仓库的 `.github/scripts/test-wheel.sh` 会复制测试到临时目录并隔离安装；CI 还会从源码包重建 wheel。平台矩阵、artifact 名称、触发条件和 PyPI 发布流程以 [CI 文档](ci.md)及 workflow 为准。

## 发布前核对

1. `pyproject.toml` 与 Cargo workspace 的版本一致，两个 lockfile 与源码同步。
2. 在实际待发布 commit 上检查 CI、wheel / sdist、隔离安装和版本元数据；release 构建模式不能替代这些结果。
3. 文档示例与公开 API、CLI 参数一致，文档链接可在干净 checkout 中解析。
4. 对外的性能结论附 fixture、环境、源码 / 二进制标识和计时协议，不把 warm-cache 合成结果推广为生产吞吐。
5. 分发前核查许可证与第三方源码说明。`refs/` 的测试依赖与真正进入扩展的依赖分别核对。

本页是执行流程，不表示这些检查已经在当前会话或当前 release 上通过。发布动作见 [CI 与 PyPI 发布](ci.md)。

## 文档与本地记录

公开文档由 [docs/README.md](README.md) 导航；`development/` 中的提议会明确标注实现状态。脚本入口与对照环境单独记录在 [scripts/README.md](../scripts/README.md)。本地优化轮次、临时日志和历史审查记录保留在被 Git 忽略的目录，不作为当前公开接口或发布证据。

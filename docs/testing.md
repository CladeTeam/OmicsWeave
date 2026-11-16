# 测试

Rust 测试按当前行为组织。格式读取、写入、查询、缓存和执行器分别拥有具名场景；错误输入与失败恢复放在对应行为旁。Python 绑定测试位于 `src/tests/`，负责公开 Python 接口、NumPy 内存布局、结果所有权和进程生命周期。

## 文件布局

小型 crate 使用少量行为命名的顶层测试文件。部分 crate 使用 `tests/integration.rs` 注册 `cases/` 中的测试；`convert`、`vcf` 等则按行为拆成多个顶层 target，不能仅运行名为 `integration` 的 target 代替整个 crate 的测试。多数 crate 的 `support/` 保存共用夹具、故障注入或独立预期模型；`quant` 将共用夹具、独立 oracle 和输出读取放在 `common/`，`support/` 负责跨线程数、spill 与存储格式的执行轴检查。需要全局分配器、C 全局状态或子进程入口的测试可以使用独立 target。

只为真实共享的逻辑建立 helper。测试的预期结果应来自明确的字面值、简单独立模型或参考实现；单独的读写往返不能证明文件格式兼容。计时和吞吐测量放在 `benches/`，不混入默认正确性测试。

## 运行

完整测试会编译 C/C++ 参考库和整个 Rust workspace，耗时较长且占用较多内存，建议在资源充足的构建机上运行。

```sh
# 在项目根目录运行。
cargo test -p omicsweave-core --locked
cargo test --workspace --locked
cargo test --workspace --release --locked
cargo fmt --all -- --check
cargo check --workspace --all-targets --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
```

首次构建参考库会编译 C/C++ 和 HDF5，构建目录应有充足空间。并行修改期间，使用独立的 `CARGO_TARGET_DIR` 和受控的 `CARGO_BUILD_JOBS`，避免干扰其他任务。依赖已经缓存时可以加 `--offline`。

Python 测试必须使用当前源码构建的扩展，或者在独立环境中使用待验证 wheel；不能将旧扩展的结果作为当前 Rust 代码的验证。

```sh
uv run maturin develop --locked
uv run ruff check src .github/scripts
uv run ruff format --check src .github/scripts
uv run pytest src/tests
```

CI 的具体工具链、Ruff 版本和执行条件以 [ci.yml](../.github/workflows/ci.yml) 为准：最低支持 Rust 版本运行 `cargo check --workspace --all-targets`，指定工具链运行 Clippy，stable 运行 Rust 测试；Python lint 和格式检查由独立 lint job 执行。

## 用例与夹具

- 每个公共行为至少检查正常路径和适用的空输入、边界、非法输入。
- 对调用者缓冲区检查精确长度、offset、字段顺序以及前后守卫；验证失败后输出和内部状态的约定。
- 对有状态 reader、writer、cache、executor 检查连续调用、clone、预算不足、错误传播和资源生命周期。
- 对并行、分块或外存排序路径，使用相同输入比较不同线程数、分块和 spill 设置的结果。
- 参考实现和被测实现必须分别执行。参考工具不可用、未定义输入和实现失败应分别报告，不能把未比较当成比较通过。
- 外部工具测试单独标明运行条件；运行该组测试前检查依赖齐全。不要使用整仓库 `--ignored` 来选择它们，因为 ignored 也可能包含需要特定父进程的入口或真实大数据检查。

## 属性测试

属性测试使用小而有代表性的结构化输入，尽量允许缩减到具体记录、区间或矩阵。失败报告包含输入和影响行为的参数。

保存的 proptest 反例属于测试输入。移动测试时同步迁移持久化配置并确认反例会重放；也可以把最小输入整理成具名的确定性场景。文件路径和测试入口变化后，不得假设旧 persistence 路径仍然有效。

## 验证结果的边界

报告应分别记录正确性测试、参考工具对照、Python 绑定和性能测量的执行情况。明确列出失败、跳过及未运行项目。覆盖率用于定位尚未检查的分支；测试数量、单纯往返成功和覆盖率百分比都不能代替行为覆盖清单。

## 模块入口

下表列出当前顶层测试 target；部分入口注册 `cases/`，另一些直接包含具名场景。共享辅助代码的职责按上述文件布局区分，尤其注意 `quant/common` 与 `quant/support`。

| crate | 测试入口 |
|---|---|
| bam | [integration](../crates/bam/tests/integration.rs) |
| bbi | [integration](../crates/bbi/tests/integration.rs), [resource_contracts](../crates/bbi/tests/resource_contracts.rs) |
| bgzf | [index_oracle](../crates/bgzf/tests/index_oracle.rs), [integration](../crates/bgzf/tests/integration.rs) |
| consensus | [integration](../crates/consensus/tests/integration.rs), [single](../crates/consensus/tests/single.rs) |
| convert | [alignments](../crates/convert/tests/alignments.rs), [contacts](../crates/convert/tests/contacts.rs), [formats](../crates/convert/tests/formats.rs), [pipeline](../crates/convert/tests/pipeline.rs) |
| cool | [integration](../crates/cool/tests/integration.rs), [real](../crates/cool/tests/real.rs) |
| core | [exec_lifecycle](../crates/core/tests/exec_lifecycle.rs), [integration](../crates/core/tests/integration.rs) |
| cram | [integration](../crates/cram/tests/integration.rs) |
| fasta | [integration](../crates/fasta/tests/integration.rs) |
| ffm | [integration](../crates/ffm/tests/integration.rs) |
| gtf | [integration](../crates/gtf/tests/integration.rs) |
| hic | [integration](../crates/hic/tests/integration.rs) |
| htscodecs-sys | [codecs](../crates/htscodecs-sys/tests/codecs.rs) |
| motif | [external_tools](../crates/motif/tests/external_tools.rs), [integration](../crates/motif/tests/integration.rs) |
| msa | [integration](../crates/msa/tests/integration.rs) |
| python | [Python 领域测试](../src/tests/) |
| quant | [alleles](../crates/quant/tests/alleles.rs), [atac](../crates/quant/tests/atac.rs), [counts](../crates/quant/tests/counts.rs), [coverage](../crates/quant/tests/coverage.rs), [feature_intervals](../crates/quant/tests/feature_intervals.rs), [hic](../crates/quant/tests/hic.rs), [labeling](../crates/quant/tests/labeling.rs), [long_indices](../crates/quant/tests/long_indices.rs), [methylation](../crates/quant/tests/methylation.rs), [splicing](../crates/quant/tests/splicing.rs), [task](../crates/quant/tests/task.rs) |
| sam | [integration](../crates/sam/tests/integration.rs) |
| vcf | [encoding](../crates/vcf/tests/encoding.rs), [fields](../crates/vcf/tests/fields.rs), [filter](../crates/vcf/tests/filter.rs), [indexed](../crates/vcf/tests/indexed.rs), [reader](../crates/vcf/tests/reader.rs) |

## 外部参考与真实数据

默认正确性测试使用随源码构建的参考库与小型确定性夹具。MEME 工具对照和真实 `.mcool` 是独立的可选入口；它们不以缺依赖后直接返回的方式显示通过。

`motif` 的四个工具比较依赖相应的 MEME 5.5.9 可执行文件。设置 `OMICSWEAVE_FIMO`、`OMICSWEAVE_MAST`、`OMICSWEAVE_MCAST`、`OMICSWEAVE_GLAM2SCAN` 后显式运行：

```sh
cargo test -p omicsweave-motif --test external_tools --locked -- --ignored
```

只验证一个工具时，按其完整测试名称选择用例，并只提供该工具的路径。显式运行缺依赖会失败。`cool` 的真实数据检查使用 `OW_COOL_SMOKE` 指定已有 `.mcool`，按 `tests/real.rs` 的条件运行；不自动获取大数据。

C 参考程序的进程退出、全局状态和不可比较输入不能影响主测试进程。VCF 与 consensus 的 helper 使用独立可执行进程，父进程验证入口握手、结果完成标记和超时；helper 自身在默认测试列表中标为 ignored，但由对应父用例实际调用。VCF 严格差分使用类型与形状明确的表达式家族，并要求参考实现完成每条记录；畸形语法和输入另行验证安全拒绝、克隆状态与错误后复用。

## 输出缓冲区

固定大小 Python 查询的 `out=` 必须满足类型、形状、连续性和可写要求。拒绝输入时不改写调用者的数组，合法数组可在失败后继续使用。原生绑定通过可失败的 NumPy 可写借用返回正常 Python 异常；只读或已被借用的数组不能触发 Rust panic。默认查询返回独立拥有的结果，继续查询或关闭 reader 不影响已有结果。

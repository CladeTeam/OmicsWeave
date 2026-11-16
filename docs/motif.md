# Motif 匹配

[FFM 格式](formats/ffm.md) · [格式转换](convert.md) · [按区域查询](query.md)

> 状态：Rust 实现位于 `crates/motif`（`omicsweave-motif`），覆盖 FIMO、MAST、GLAM2Scan 与 MCAST；Python 接口为 `omicsweave.motif`（见[Python 接口](#python-接口)）；CLI 尚未实现。

motif 匹配层把 [FFM](formats/ffm.md) 中的 motif 模型编译成扫描程序，对序列扫描一次，同时更新所有 motif 的规约状态。MEME Suite 里「一个二进制一个功能」的几个匹配工具，在这里是同一套组件的不同组合：

| 工具组合 | 模型 | 打分器 | 规约 | 对应文本接口 |
|---|---|---|---|---|
| FIMO | 因子模型 | 分块查表 | 过阈值的窗口 | `fimo::write_text` |
| MAST | 因子模型，预设面向 PWM | 分块查表 | 最佳窗口、非重叠排布、序列排名 | `mast::write_hit_list` |
| MCAST | 因子模型，预设面向 PWM | 分块查表 | 有序命中聚簇、尾部统计 | `mcast::write_text` |
| GLAM2Scan | 带 gap 的 profile | 局部比对 DP | 比对流、每个 profile 的全局前 N | `profile::write_glam2_text` |

这些预设以 MEME Suite 5.5.9 为对照目标；可选黑盒测试及其运行条件见[验证与性能证据](#验证与性能证据)。实现了预设、测试存在或测试提前返回，都不等于当前环境完成了外部工具对照。

## 一个完整的 Python 例子

下面的两行计数分别偏好 A 和 C，模型窗口长 2；输出坐标始终相对于输入序列、从 0 开始且右端不含。

```python
from omicsweave import motif, types

dna = types.Alphabet.dna()
model = types.Motif.from_pfm("AC", [[9, 1, 1, 1], [1, 9, 1, 1]])
program = motif.compile(
    [model], motif.Background.uniform(dna), motif.Scoring.exact(),
    alphabet=dna, strands="+", threads=1,
)
hits, counts = program.scan(
    [b"AACGT", b"CCAA"], motif.Hits(score=1.0), motif.Counts(score=1.0), threads=1,
)
# hits 的每一列长度相同；counts 是 2 条序列 × 1 个 motif 的 CSR 计数。
print(hits.start, hits.end, counts.shape)
```

传入 `.ffm` 路径或 `query.MotifReader` 时从文件取得字母表；传入 Python `Motif` 值时必须显式给出 `alphabet=`。编译持有独立的打分表，之后读取器关闭不会使程序失效。

## 分层

```text
工具层   fimo · mast · mcast · profile（GLAM2Scan）       只做组装：打分预设 + 规约 + 统计 + 文本输出
统计     stats：BestOf · Combine · Adjust · TailFit         可替换的接口；MEME 的公式各是一种实现
规约     窗口级（无序、线程本地）→ 序列级（按序列、按位置有序）→ 集合级（finish）
打分器   分块查表（因子模型）· 局部比对 DP（带 gap 的 profile）
模型     FFM 因子模型 / 带 gap 的 profile；Background；Scoring（含列先验、模糊字母规则）
```

## 责任边界

核心（模型、打分器、规约框架）只负责**一个窗口在声明的背景模型下的精确统计量**，以及对这些量做确定性的汇总。把窗口统计量变成序列级或集合级统计量的每一步——最佳窗口的序列 p 值、多个 motif 的组合、E 值与 q 值、对抽样分数的尾部拟合——都是 `stats` 中的接口，调用方可以换成自己的实现（闭包也可以）；MEME Suite 的算法只是各接口的一种实现，供工具层复现 MEME 的输出。

| 框架负责 | 框架不负责 |
|---|---|
| 背景模型：显式构造、截断阶数、对称化、从序列估计、生成合成序列 | motif 文件格式解析（`ow convert` → `.ffm`） |
| 打分：伪计数或列先验、log-odds、模糊字母、正反链 | 读取序列文件、解析 FASTA 头（调用方给出内存中的序列、名称与坐标原点） |
| 量化、零分布、单窗口 p 值、阈值 ↔ 分数换算 | 选定某一种多重检验或富集统计作为唯一做法 |
| 编码、扫描、剪枝、并行、序列级有序交付 | motif 发现、motif 比较、结果的生物学解释 |
| 确定性规约与各工具的组装、MEME 兼容的文本输出 | |

## 总体流程

```text
FFM MotifView ─┐
Background ────┼─► ① 估计 ─► ② 定向 ─► ③ 对比 ─► ④ 量化 ─┬─► ⑤ 零分布 / 阈值 ─┐
Scoring ───────┘   概率因子   每条链    ScoreModel   IntModel └─► ⑥ 降为内核 ─────┼─► ⑦ 链接 ─► 扫描
规约（Reducer）的需求 ──────────────────────────────────────────────────────────┘
```

| 步骤 | 输入 → 输出 | 公开产物 |
|---|---|---|
| ① 估计 | FFM 频率 → 归一化、平滑后的全局编码概率 | — |
| ② 定向 | 每条链一份 | — |
| ③ 对比 | 加上背景项 → 加性项 | `ScoreModel` |
| ④ 量化 | 整数项 + 还原参数 | `IntModel` |
| ⑤ 零分布 | 整数分数的精确分布；p 值 ↔ 分数 | `NullDist`、`Variant` |
| ⑥ 降为内核 | 分块查表、剪枝顺序 | — |
| ⑦ 链接 | 规约需求 → 每个（motif，链）的内核模式 | — |

①–⑥ 在 `Program::compile` 中完成，⑦ 在每次 `Program::scan` 开始时完成。单个 motif 的 ①–⑤ 也可以用 `MotifModel::compile` 单独完成，用于阈值标定。

## 输入

### 模型

FFM 的 motif（`MotifView`）：宽度 `W`、support、可选 `sites`、因子列表。匹配层只读，不修改格式。

### 背景

背景必须显式给出，没有默认值。`Background` 是 `m` 阶 Markov 模型（`m ≥ 0`），由 `k = 1 … m+1` 的正 k-mer 频率定义（按给定值保存，无需归一化）；`c ≤ m` 个上下文符号之后的条件概率是 (c+1)-mer 频率对最后一个符号归一化的结果。

| 构造 | 说明 |
|---|---|
| `uniform(alphabet)` | 0 阶均匀 |
| `nrdb(alphabet)` | MEME 内置的 NRDB 0 阶频率（`ACGT` 与 20 种氨基酸） |
| `frequencies(p)` | 0 阶频率 |
| `from_kmer_frequencies(n, kmers)` | 1…m+1 阶 k-mer 频率 |
| `from_meme(alphabet, text)` | MEME 背景文件的文本（`fasta-get-markov` 输出） |
| `estimate(alphabet, sequences, m, pseudocount)` | 从内存中的序列计数估计 |

变换：Rust `symmetrized(complement)` 把每个 k-mer 频率与其反向互补平均；`truncated(c)` 截断到 c 阶。Python 对应 `bg.symmetrized(alphabet)`、`bg.truncated(order)`，MEME 文本入口为 `Background.parse_meme(text, alphabet)`。频率要求有限且严格为正，每阶数组分别有 `n^k` 个元素，至多 `2^24` 格；`estimate` 与 `Program` 所用的准备后背景阶数至多 32，表规模可能使实际可用阶数更低；Python 值对象的构造不替代原生编译时的完整校验。

Rust `sample(&alphabet, &lengths, seed)` 返回 `(Vec<i64>, Vec<u8>)`；Python `bg.sample(alphabet, lengths, seed=0)` 返回 `StringVector`。每条序列重新开始 Markov 上下文；长度可以为 0。Python 长度必须是一维非负整数，背景与字母表的符号数必须一致，总长度须适合地址空间与 `i64` 偏移。Rust 的非 `Result` 采样接口对后两项违反前提的情况提前 panic；Python 包装层先校验输入。采样按 seed 确定，使用滚动上下文，不保留整条序列的整数编码副本。

### 打分策略 `Scoring`

下表的默认指 Rust `Scoring::default()`，即 `exact()`；**Python `motif.compile(..., scoring=None)` 默认采用 `Scoring.fimo()`**，两者不可混淆。

| 字段 | 含义 | `Scoring::fimo()` | `Scoring::exact()`（默认） |
|---|---|---|---|
| `pseudocount` | 先验总质量 α | 0.1 | 0.1 |
| `default_sites` | motif 没有 `sites` 时的观测数 | 20 | 20 |
| `truncate_sites` | `sites` 取整数部分 | 是 | 否 |
| `normalize_tolerance` | 和与 1 相差不超过它的归一化组原样使用 | 1e-5 | 0 |
| `log` | 对数 | `MemeBits`：`ln(x) · 1.44269504` | `Bits`：`log2(x)` |
| `quantizer` | 见 ④ | `Fimo { range: 1000 }` | `Precision { states: 2^16 }` |
| `rounding` | 报告的分数与 p 值的舍入 | 10 位有效数字（MEME `RND`，就近偶数） | 无 |
| `background` | 背景的准备 | `Meme` | `AsGiven` |
| `null_budget` | 单个零分布允许的最大计算步数 | 2^33 | 2^33 |
| `null_letters` | 零分布专用的字母概率（原样使用，不归一化） | 无 | 无 |
| `shared_reverse_null` | 反链沿用正链的零分布（每个 motif 一张 p 值表） | 否 | 否 |
| `column_prior` | 代替伪计数的单位置列先验（`ColumnPrior`） | 无 | 无 |
| `ambiguity` | 模糊字母及其打分规则（`Ambiguity`） | 不打分 | 不打分 |

`log` 还有 `Natural` 与 `Centibits`（`ln(x) / ln 2 · 100`，MAST 的单位）；`quantizer` 还有 `Mast { range }`（offset 取最小值本身、scale 为 `range / (max − min)` 保留 14 位有效数字）。

**模糊字母**（`Ambiguity`）：字母表之外、代表一组符号的字节（IUPAC 的 `R Y S W K M B D H V N` 等）。规则 `Skip` 时含它们的窗口不打分（FIMO）；`Mean` 时取成员符号分数按零分布字母概率加权的平均（MAST）；`Sum` 时取和（MCAST）。只有单位置项为模糊字母建格，格子参与量化的取值范围；含模糊字母的窗口用同一张（只含字母表符号的）零分布取 p 值。

**列先验**（`ColumnPrior`）：一个 trait，把单位置 joint 因子的归一化频率映射为平滑后的概率，代替伪计数，例如 MCAST 的替换矩阵先验。输入使用局部编码；有 rest 时最后一个元素为 rest，长度为局部 arity。返回表必须长度相同、值有限且非负、至少一个值为正，否则编译返回错误。先验输出不被再次归一化，调用方负责它的概率语义。

Python 对应 `Scoring.column_prior`：`SubstitutionPrior` 或 `callback(column, sites, letters)`。回调的 `column`、`letters` 为一维 `float64` 数组，返回同长的一维概率数组；通过 `Scoring.replace(column_prior=...)` 配置。调用自定义 Python 逻辑时会重新取得 GIL。

`BackgroundRule::Meme` 按 MEME 的方式准备背景：截断为 0 阶；扫描双链时与互补平均；每个字母加 5e-7 后除以总和，再调用容差为 0 的 `meme_normalize`。准备后的字母频率由条件概率表原样保存。

### 链

`Strands::Forward` / `Reverse` / `Both`；字母表没有互补表时只能 `Forward`。

### 序列

序列以内存中的 ASCII 字节给出，例如 `SequenceReader.fetch(...).bases` 中的碱基：Rust 单条用 `Program::scan_sequence(&[u8])`，批量用 `Strings`（Arrow `LargeString` 布局的 `i64` 偏移 + 数据），每行一条序列。字节经字母表的 LUT 映射为全局编码；LUT 为 `0xFF` 的字节是未知符号。序列名称、坐标原点等由调用方保存，只在输出层使用。

## 打分语义

记窗口起点为 `i`、motif 位置为 `0 … W-1`，窗口内符号为 `x = s[i … i+W)`。正链分数是对数似然比：

```text
S(i) = M(x) − B(x)

M(x) = Σ_f  m_f(x_{P_f})                        motif 部分，按因子求和
B(x) = Σ_{p ∈ C}  ln b(s[i+p] | s[i+p−c … i+p))  背景部分，每个覆盖位置一次
```

- `m_f`：joint 为 `ln q̃(x_P)`，conditional 为 `ln q̃(x_t | x_ctx)`，potential 为 `ln v`。`q̃` 是平滑后的频率（见「伪计数」）。
- `C`：模型生成的位置——joint 因子的全部位置、conditional 因子的最后一个位置。只被 potential 覆盖或完全不被覆盖的位置不贡献，等价于「motif 在这些位置与背景相同」。
- **背景上下文取真实序列。** 位置 `i+p` 的条件概率读取它之前 `c = min(m, 可用长度)` 个符号，可以越过窗口左端。可用长度在序列开头或未知符号处截断，此时使用 c 阶条件概率，等价于 Markov 链在该处重新开始。
- **窗口有效性与 FIMO 相同。** 窗口内任一符号未知则不打分；上下文中的未知符号只缩短上下文。
- **反链。** 反链分数为 `M(rc(x)) − B(x)`：motif 部分取反向互补，背景部分仍按正链读取，因为被检验的是同一段 DNA 在背景下的概率。FIMO 的做法相同（用正链背景去除反链 motif 的频率）。
- **因子语义。** PWM 与 conditional 链（BaMM）的 `M` 是窗口的真实对数概率，`S` 是严格的似然比。重叠的 joint 因子（二核苷酸模型）不构成归一化分布，`M` 是能量式的打分；p 值由该分数的精确分布给出，因此仍然精确。

PWM、0 阶背景时 `S = Σ_p log(q̃_p(x_p) / b(x_p))`，即 FIMO 的 PSSM。

## 伪计数

记 α 为 `pseudocount`，`n` 为 `sites`（缺失时为 `default_sites`，`truncate_sites` 时取整数部分），`b₀` 为背景的字母分布。`p̂` 为按归一化组归一化后的频率。

| 因子 | 平滑 | 先验均值 |
|---|---|---|
| joint，1 阶 | `q̃ = (n·p̂ + α·b₀) / (n + α)` | `b₀` |
| joint，k ≥ 2 阶 | `q̃ = (n·p̂ + α·μ) / (n + α)` | `μ(x_P) = Π_{p∈P} q̃_p(x_p)`，`q̃_p` 为该因子在位置 p 上的边缘按 1 阶规则平滑 |
| conditional | 每个上下文：`q̃(·|ctx) = (n·p̂(·|ctx) + α·b₀) / (n + α)` | `b₀` |
| potential | 不平滑 | — |

理由：

- 1 阶规则即 FIMO（MEME `apply_pseudocount_to_motif`），先验均值是背景。
- 高阶 joint 向「自身边缘的乘积」收缩，而不是向背景收缩：观测不足时退化为位置独立的 PWM，只平滑位置间的耦合。若向背景收缩，两个各自高频、但组合未被观测到的符号会被重罚。
- conditional 表不携带上下文的出现频率，无法得到每个上下文的真实证据量；每个上下文按 `n` 计，作用是保证概率为正、尽量少改动来源已平滑好的概率（BaMM 等来源自带层次平滑）。
- potential 是能量，0 表示禁止，不做平滑。

**support 与 rest**：平滑在局部编码上进行，rest 的先验质量为 `b₀(rest) = Σ_{c∈rest} b₀(c)`。展开到全局编码时，rest 中符号 `c` 在每个生成位置的概率为 `q̃(rest) · b₀(c) / b₀(rest)`，即 motif 在 rest 内部与背景同分布。0 阶背景下对数似然比为 `ln q̃(rest) − ln b₀(rest)`，与格式文档的参考语义一致。

**α = 0**：不平滑；频率为 0 的组合记为禁止，窗口命中禁止组合时没有分数，零分布中单独计入「禁止」质量。含禁止组合的输入不应仅凭预设名称推断与外部工具完全等价。

## 编译

### ① 估计

对每个因子归一化（和与 1 相差超过 `normalize_tolerance` 时除以和）、平滑，并展开到全局编码。有限大计数的和溢出时按组内最大值缩放；观测数与伪计数的和溢出时按权重比例平滑。普通数值路径保留预设的算术顺序。

### ② 定向

反链：位置 `p → W−1−p`，每个轴按互补表重排，再把位置重排为递增顺序。

### ③ 对比

对每个覆盖位置 `p` 加入背景项。上下文长度为 0（0 阶背景，或截断到 0 的上下文）时，若 `p` 上有单位置的 joint / conditional 因子，背景并入它，单元格为 `log(q̃ / b)`，保留预设对比值整体取对数的顺序；否则背景项单独成项，位置为 `{p−c, …, p}`，表为 `−log b(x_p | 前 c 个)`。

规范化：位置集合相同的项合并；位置集合包含于另一项的项并入后者；按最后一个位置排序。结果是 `S = Σ_t T_t[x_{P_t}]`，最小位置 `lo ≤ 0` 表示读取窗口左侧的上下文。

每条链为每个上下文长度 `c = 0 … m` 各生成一份模型（`StrandModel::variants()[c]`），`c = m` 是普通窗口所用的模型，其余用于序列开头或未知符号之后的窗口。

### ④ 量化

整数打分让扫描与执行顺序无关，也让零分布是精确的离散分布。p 值对整数分数精确。

**`Fimo { range }`**（复刻 MEME `scale_score_matrix`，每条链各算一次）：

```text
small, large = 所有项所有有限格子的最小、最大值；相等时 small = large − 1
offset = floor(small)
scale  = floor(range / (large − offset))
cell   = nint((x − offset) · scale)               nint：0.5 远离零
raw(S) = S / scale + T · offset                   T 为项数，PWM 时 T = W
```

**`Precision { states }`**：每个项用自己的最小值作 offset，使所有格子非负且范围最紧；`scale = (states − 1) / Σ_t (max_t − min_t)`；`raw(S) = S / scale + Σ_t offset_t`。舍入误差 `|raw(S) − S_real| ≤ T / (2·scale)`。

截断上下文的模型沿用同一条链完整模型的 `scale`，各自取 offset，分数在实数尺度上可比。禁止格子记为 `FORBIDDEN`。

### ⑤ 零分布与阈值

零分布是窗口（连同它读取的上下文）在背景下随机生成时整数分数的分布，Markov 链从读取范围的第一个位置开始，以可用的最长上下文生成，与截断上下文窗口的打分语义一致。

- **单位置项 + 0 阶背景**：逐位置卷积，加法顺序与 MEME `get_pdf_table` 相同（外层位置、中层字母、内层旧分数升序）；生存函数从高分向低分累加，每步与 1 取最小值。正反链各算一次；外部 p 值表的一致性由对应黑盒用例检查。
- **一般情形**：沿读取范围做变量消元。位置 `p` 的状态是之前仍会被读取的位置上的符号——被结束于 `p` 及之后的项引用，或作为背景上下文。相邻位置的高阶模型状态很小；远距离的位置对只多一个符号，不会把整个跨度带入状态。
- **一般变量消元的预算**：估算步数超过 `null_budget`，或状态内存超过 2^27 个格子时编译报错，不静默改用近似。0 阶单位置卷积使用独立路径，当前不通过该步数预算门控。

阈值换算：`cutoff(p) = min { S : rounding(sf[S]) ≤ p }`，`sf` 单调，二分查找。没有分数能达到阈值时，该（motif，链）不产生普通窗口；扩展的模糊字母分数仍单独判断，可能高于普通分数的上界。`Variant` 还给出 `p_min` 与 `alpha_eff(threshold)`，可直接用于标定。

### ⑥ 降为内核

**编码约定**：每个符号占 `bits = ceil(log2 |A|)` 位（DNA 2、蛋白 5），k-mer 编码为 `k` 个符号依次拼接，`k` 取满足 `2^(bits·k) ≤ 1024` 的最大值（DNA k = 5，蛋白 k = 2）。

**分块**：块是读取窗口偏移 `s` 处 k-mer 前缀的一张表，表值为分给它的项之和。位置跨度 `[a, b]`（`b − a < k`）的项能被起点 `s ∈ [b−k+1, a] ∩ [lo, W−k]` 的块读到；按可行区间右端贪心穿刺，块数最少。`W − lo ≤ k` 时只有一个短块。跨度超过 `k` 的项成为稀疏块，直接按 `code[i+p_j]` 计算下标。

| 模型 | 块 |
|---|---|
| PWM，W = 12，DNA | 偏移 0、5、7（最后一块右对齐） |
| 相邻二核苷酸 | 块起点以 k − 1 为步长 |
| 1 阶 Markov 背景下的 PWM | 背景项跨两个位置，偏移从 −1 开始 |

含禁止格子的单元格记为 `−(max + 1)`，任何含禁止组合的窗口总分为负。

**剪枝**：块按表值的跨度降序排列，预计算剩余块的最大和 `rest_max[j]`。前 `split` 块对整个 tile 稠密求和；部分和加 `rest_max[split]` 达不到阈值的窗口直接跳过，其余逐块补全并在上界不足时提前停止。`split` 在每次 scan 的链接阶段根据背景随机窗口的存活率估计代价后选择，只影响速度，不影响结果。

## 执行

- **链接**：每个规约声明需求（`Need::All` / `Above(threshold)` / `Best` / `Ordered(threshold)`），可以用 `Demand { need, part }` 限定到程序的某个部分；对每个（motif，链）取并：有 `All` 即为稠密模式；否则为阈值模式，阈值取各 `Above` 与 `Ordered` 的最小整数 cutoff；`Best` 另外跟踪每条序列的最佳窗口，以当前最佳作为动态阈值。没有需求的（motif，链）不扫描。
- **多个部分**：`Program::combine` 把多个 `Program`（各自的背景、打分策略与链）合成一个，motif 编号依次相接，序列只编码一次。规约用 `in_part(i)` 只订阅其中一部分。
- **tile**：序列装入 tile（默认约 64 K 符号）。长序列切成多段，每段带上最多 `m` 个上下文符号和最后一个窗口读取的碱基；短序列拼进同一个 tile，以一个未知符号分隔，跨序列的窗口因含未知符号自然无效。
- **编码一次**：每个 tile 计算一次全局编码、滚动 k-mer、从每个位置起的已知符号连续长度、之前的已知符号数和窗口归属，所有 motif 共用。
- **边缘窗口**：上下文不足 `m` 的窗口（每条序列开头及每段未知符号之后至多 `m` 个）在编码时列出，用对应上下文长度的模型逐项求和。
- **模糊字母窗口**：打分策略为模糊字母建格时（仅 0 阶背景），含模糊字母、不含未知字节的窗口在编码时按模糊字母的位置列出，逐项求和。
- **有序阶段**：有 `Ordered` 需求时，各 tile 中通过阈值的窗口送到一个协调线程，按 tile 顺序重新排好，再按序列、按起点交给 `Reducer::sequence`；长序列分多段依次交付，状态可以跨段携带。工作线程最多领先协调线程固定数量的 tile，缓冲有上界。
- **并行**：工作单元是（tile，motif 组），输入很小时也能用满线程；规约在线程本地累积，结束时合并。`threads = None` 为可用并行度。
- **确定性**：结果与线程数、tile 大小无关。

## 规约

```rust
pub trait Reducer: Send + Sized {
    type Output;
    fn needs(&self, out: &mut Vec<Demand>);
    fn prepare(&mut self, program: &Program, info: &ScanInfo<'_>) -> Result<()>;
    fn fork(&self) -> Self;
    fn windows(&mut self, batch: &Windows<'_>);                         // 窗口级：线程本地、无序
    fn sequence(&mut self, program: &Program, chunk: &SequenceChunk<'_>); // 序列级：按序列、按位置有序
    fn merge(&mut self, other: Self);
    fn finish(self, program: &Program, info: &ScanInfo<'_>) -> Self::Output; // 集合级
}
```

三级规约：

- **窗口级**：一批窗口给出整数分数与上下文长度，`Windows` 负责换算报告分数、p 值与阈值判断；同一序列的窗口可能分散在多批中，顺序不定。
- **序列级**：`SequenceChunk` 是一条序列的一段（`first` / `last` 标出首尾），其中通过 `Ordered` 阈值的匹配（`Match`：分数、p 值、motif、起止、链）按起点、motif、链排序，跨所有 motif 与两条链。长序列可以多次交付，但各段按顺序覆盖起点范围，每个起点只归属一段；空序列没有有序段回调。
- **集合级**：`finish` 拿到全部输入（序列数与长度），用于 E 值、排名、尾部拟合等。

元组 `(A, B)` … `(A, B, C, D)` 本身也是规约，一次扫描驱动多个规约。

| 规约 | 需求 | 输出 |
|---|---|---|
| `Hits::new(threshold)` | `Above` | `Vec<HitEntry>`：分数、p 值、序列、motif、起止（0-based 半开）、链；按 motif、序列、起点、链排序 |
| `SequenceBest` | `Best` | `Vec<BestEntry>`：每个（序列，motif）在双链上的最佳窗口；平局取起点更前、正链优先 |
| `SequenceCounts::new(threshold)` | `Above` | `CountMatrix`：序列 × motif 的稀疏计数（CSR），双链合计 |
| `Dense` | `All` | `DenseScores`：每条序列 `[motif][链][位置]` 的报告分数，无分数处为 NaN |

`HitEntry` 与 `BestEntry` 是 `#[repr(C)]` 行，64 位构建中分别为 40 与 32 字节、自然对齐。Python 扫描接口返回各字段的独立 NumPy 列，不要求调用方构造这两种结构化 dtype。

### 公开模型与编译边界

`ScoreModel`、`Term`、`IntTerm` 与 `IntModel` 的字段公开，用于检查或独立计算；它们不是可写的 `Program` 视图。修改导出的副本不修改已编译程序。公开 `IntModel::score` / `score_ambiguous` 使用安全的索引和整数运算：位置不足、未知编码、缺少表格值、空位置项或累计分数超出 `i32` 时返回 `None`。它们不因为模型来自某次编译就信任调用方后续修改。

`Program` 内部的查表内核仅处理编译时校验并私有保存的布局，局部 `unsafe` 的边界依据写在 `SAFETY` 注释中。频率解析、元信息和上述公开可修改模型的校验保留 safe 实现。编译还要求 `width + background.order() ≤ i32::MAX`，全局展开表最多 `2^24` 格；合法 FFM 文件也可能因匹配层的这些限制而不能编译。

`NullDist` / `NullDistribution` 的 `pmf` 只包含有有限分数的质量，`forbidden` 另存。`pvalue(s ≤ 0)` 返回 `survival[0]`，所以存在禁止质量时可以小于 1；大于普通分数上界时返回 0。

## FIMO 兼容

`Scoring::fimo()`、`Hits::new(Threshold::PValue(p))` 与 `fimo::write_text` 组合起来，面向以下输出配置：

```sh
fimo --text --no-qvalue --bfile <背景> --thresh <p> motifs.meme seqs.fa
```

（`--norc` 对应 `Strands::Forward`；`--bfile --uniform--` / `--nrdb--` 对应 `Background::uniform` / `nrdb`。）

| 项目 | 处理 |
|---|---|
| 背景 | 截断为 0 阶；双链时与互补平均；加 5e-7 后归一化 |
| motif 行 | 和与 1 相差超过 1e-5 才重新归一化 |
| `sites` | 取整数部分，缺失为 20 |
| 对数 | `ln(x) · 1.44269504`，对 `q / b` 整体取对数 |
| 打分表与 p 值表 | 正反链各自量化、各自卷积 |
| 判定 | `RND10(sf[S]) ≤ thresh` |
| 分数输出 | `RND10(S/scale + W·offset)`，`%g` |
| p 值输出 | `RND10(sf[S])`，`%.3g`；q 值列为空 |
| 坐标 | 1-based 闭区间，加上调用方给出的坐标原点 |
| 匹配序列 | 正链原样输出；反链输出字母表的互补符号（规范大写） |
| 行顺序 | motif → 序列 → 起点 → 正链、反链；没有命中时不输出表头 |

不复刻的行为：FASTA 头解析（`--parse-genomic-coord` 等，名称与坐标原点由调用方给出；FIMO 会把字母表外的字符替换为 `N`，调用方需要同样处理时自行替换）；`--max-strand`；非 `--text` 模式下 `--max-stored-scores` 的丢弃与基于抽样的 q 值；`--psp` 先验；α = 0 时与外部工具对零概率的特殊处理。

## MAST

`mast::compile`（`mast::scoring` 预设）与 `mast::Mast` 规约组合起来，对应 `mast`；`mast::write_hit_list` 输出 `-hit_list` 的行。

| 组件 | 处理 |
|---|---|
| 打分 | 伪计数 0.01、缺省 100 个 sites（取整数部分）；`log(f/b)/log(2)·100`；`Quantizer::Mast { range: 100 }`；打分背景按 MEME 方式准备（MAST 默认用 motif 文件的背景，由调用方给出） |
| p 值 | 每个 motif 一张表（反链共用），字母概率为 MAST 内置频率或背景文件原值，保留 8 位有效数字、双链时与互补平均，不归一化 |
| 模糊字母 | IUPAC（另有 `X`、`.` 代表 `N`），取成员分数按 p 值字母概率加权的平均 |
| 每个 motif 的最佳窗口 | `Need::Best`；按可还原的 raw score 比较，平局取起点更前、再正链优先，保留该窗口自身 context/strand 的 p 值；`Combine` 窗口数加倍，`Separate` 分别计，`Forward` 只看正链 |
| 序列 p 值、组合、E 值 | `MastOptions` 中的 `best_of`（默认 `p·n` / `1 − (1 − p)^n`）、`combine`（默认乘积的 p 值）、`adjust`（默认 p × 序列数） |
| 排布 | `Need::Ordered`：逐 motif、逐起点贪心，p 值小于所覆盖的已选匹配 p 值之积时替换它们 |
| 排名 | 先按 MAST 的预筛（p × 已读序列数 < 阈值），再按 E ≤ 阈值；按 p 值、输入顺序、链排序 |
| 相关系数 | `mast::correlations`（整数打分列的 Pearson 相关，取最佳重叠的均值）与 `mast::redundant`（`-remcorr` 规则） |


不复刻的行为：motif 文件中自带的 `log-odds matrix` 段（FFM 只存频率，分数总是由频率计算）；`-comp`；超过 100,000 个碱基的序列在 `-hit_list` 中的分段排布；`-diag` 的间距项；序列中出现非法字符时 MAST 停止读取整个库。

## MCAST

`mcast::scoring` 预设、`mcast::Clusters` 规约与 `mcast::search` 组合起来对应 `mcast`；`mcast::write_text` 输出 `--text` 的 TSV 行。

| 组件 | 处理 |
|---|---|
| 打分 | 不加伪计数；列先验为替换矩阵先验（`substitution::SubstitutionPrior`：PAM-1 转换 / 颠换矩阵、Karlin λ、`(20·f + 1·g) / 21`），实现 `ColumnPrior`；MEME bits；`Fimo { range: 1000 }`；MEME 背景准备（MCAST 默认用 NRDB，由调用方给出） |
| 模糊字母 | IUPAC，取成员分数之和，参与缩放范围；其余任何未知字节按 `N`（`Ambiguity::with_unknown`） |
| 命中 | `Need::Ordered(PValue(pthresh))`；命中分数 `−log2(p / pthresh)` |
| 聚簇 | 沿位置的动态规划：gap 代价 `D / max_gap` 每位置，簇的阈值 `D` 由 `cluster_threshold` 按 MCAST 公式给出（可直接指定）；两条链为各自独立的 motif 状态；长序列跨段携带状态，MCAST 的分段规则一并复现 |
| 统计 | `Statistics { null, tail, evalue, qvalue }`：零分布来源 `ClusterNull`（`MemeSynthetic` 复现 MCAST 的 MT19937 合成序列；`BackgroundNull` 用 `Background::sample`），尾部拟合 `ClusterTail`（`GcBinnedExponential` 为 MCAST 的按 GC 分箱指数拟合；`Unconditional` 包装任意 `TailFit`），E 值用 `stats::EValue`，q 值用任意 `Adjust`（`MemeQValue` 为 MCAST 的做法） |
| 输出 | `ScoredCluster`（坐标、分数、命中、GC、p / E / q），`write_text` 使用 MCAST 风格的 TSV 字段 |

分段由内存中的序列长度决定，不解析 FASTA 换行。本实现保留搜得的簇；分数、p 值、GC 完全相同时按发现顺序处理。外部工具的大库截断策略与全部 CLI 选项不属于此接口的承诺。

## GLAM2Scan

带 gap 的 profile（FFM 的 gap 扩展）由 `profile::ProfileProgram` 以局部比对打分，与因子模型共用 `Background`、`Strands`、内存序列输入与并行方式。

| 组件 | 处理 |
|---|---|
| 打分 | `ProfileScoring`：删除 / 插入的 beta 伪计数（GLAM2 的 D=0.1、E=2、I=0.02、J=1）与 Dirichlet 混合先验（DNA 单组分 α=0.4；蛋白质 recode3 20 组分；也可从 `-d` 文件读入）；分数为自然对数 |
| 比对 | 每条链一次局部比对 DP；回溯取得一个匹配后屏蔽已用的（列，位置）并局部重算，产生分数不增的匹配流 |
| 规约 | `AlignmentReducer`：每条序列的各（motif，链）匹配流按序列顺序交付；`TopAlignments::glam2(n, names)` 复刻 GLAM2 的共享堆（同分按名称排序），`names: Strings` 为必填的序列名称，每个 motif 一个堆；`AllAlignments` 保留全部 |
| 输出 | `write_glam2_header` 与 `write_glam2_text`，用于 `glam2scan -t` 风格的文本 |

每条链的匹配流截到 `n` 个是精确的：同一条链上的分数不增，第 `n + 1` 个匹配不可能进入前 `n`。限制：与 GLAM2 相同，每条链的比对矩阵整体驻留内存（约 `(W+1)(L+1)·8` 字节）。不复刻：GLAM2 的自定义字母表文件；蛋白质的双链扫描（字母表没有互补表时本实现报错）。

## 接口

```rust
use omicsweave_motif::{Background, Hits, Program, ScanOptions, Scoring, SequenceBest, Strands};

let views: Vec<_> = motif_refs.iter().map(|m| m.view()).collect();
let program = Program::compile(
    &views, &alphabet, &Background::nrdb(&alphabet)?, &Scoring::fimo(), Strands::Both, None,
)?;
let (hits, best) = program.scan(
    sequences,                                  // Strings：内存中的 ASCII 序列
    (Hits::pvalue(1e-4), SequenceBest::new()),  // 一次扫描驱动两个规约
    &ScanOptions::default(),
)?;

// 多个部分一次扫描：FIMO 命中与 MAST 排名
let both = Program::combine([fimo_program, mast_program])?;
let (hits, report) = both.scan(
    sequences,
    (Hits::pvalue(1e-4).in_part(0), Mast::new(MastOptions::new(StrandMode::Combine)).in_part(1)),
    &ScanOptions::default(),
)?;

// 只做标定
let model = MotifModel::compile(&view, &alphabet, &background, &Scoring::fimo(), Strands::Both)?;
let full = model.strand(Strand::Forward).unwrap().full();
let cutoff = full.cutoff(Threshold::PValue(1e-4));
let alpha = full.alpha_eff(Threshold::PValue(1e-4));
```

## Python 接口

`omicsweave.motif` 与 `convert`、`query`、`quant` 并列，是 Rust 匹配层的组合式封装。FFM 本身按项目的分层放置：

| 层 | 名称 | 作用 |
|---|---|---|
| `omicsweave.types` | `Alphabet`、`Motif`、`Factor`、`Gaps` | 值：字母表与 motif 模型（含 ID 与按 schema 类型化的注释） |
| `omicsweave.formats` | `Ffm`、`FfmWriter` | 文件：`.ffm`，写入器 |
| `omicsweave.query` | `MotifReader`（`ow.open("x.ffm")`） | 读取：整体载入并校验，按序号或 ID 取 `Motif` |
| `omicsweave.motif` | 见下 | 匹配 |

`omicsweave.motif` 由下到上分四层，每层都可以单独使用：

1. **组件**：`Background`（频率、MEME 背景文件、NRDB / GLAM2 预设、从序列估计支持范围内的 Markov 阶数；对称化、截断、按模型采样）；`Scoring`（预设 `fimo()`、`exact()`、`mast()`、`mcast()`，字段可用 `replace` 修改）及 `Quantizer`、`Ambiguity`、`SubstitutionPrior`。列先验也可以直接传 Python 函数。
2. **编译**：`compile(motifs, background, scoring, strands=...)` 得到 `Program`；`combine(p1, p2, ...)` 把不同策略的程序合成一次扫描。`Program.model(i, strand)` 打开编译结果，得到 `CompiledMotif`：
   - 实数 log-odds 项表与整数项表；
   - `scale` / `offset`；
   - 精确零分布 `NullDistribution`；
   - `pvalue(s)`、`cutoff(p)`、`score_threshold(p)`。
   这一层可以在扫描之外自建阈值与统计。
3. **扫描与规约**：`Program.scan(sequences, *reductions)` 一次扫描返回各规约的结果（numpy 列）。可用的规约：
   - `Hits` → `HitTable`
   - `Best` → `BestTable`
   - `Counts` → `CountMatrix`（CSR）
   - `Scores` → `ScoreTracks`
   - `Mast` → `MastReport`
   - `Clusters` → `ClusterTable`
   
   规约用 `part=` 选择合并程序中的一部分。带 gap 的 profile 用 `compile_profiles` → `ProfileProgram.scan` → `ProfileAlignments`。
4. **工具**：`fimo`、`mast`、`mcast`、`glam2scan` 是上面三层的固定组合（预设打分 + 一个规约），以 MEME Suite 5.5.9 为对照目标；`motif.text` 按各工具的文本格式输出（`fimo`、`mast_hit_list`、`mcast`、`glam2scan`）。

统计在 `motif.stats` 中，每一步都可替换：

- **`BestOf`**：`Sidak`、`MemeBestOf`。
- **`Combine`**：`Product`、`MinBonferroni`。
- **`Adjust`**：`EValue`、`BenjaminiHochberg`、`MemeQValue`。
- **`TailFit`**：`EmpiricalTail`。
- **MCAST 的 `ClusterStatistics`**：零样本来源为 `MemeSynthetic` / `BackgroundNull`，尾部拟合为 `GcBinnedExponential` / `Unconditional`，另含 E 值与 q 值的校正步骤。

这些对象在原生层执行，也可以直接作用于 numpy 数组，例如 `BenjaminiHochberg()(hits.pvalue)`。形状相同的 Python 函数可以放在任何一步，原生层会回调它；回调抛出的第一个异常在调用返回时原样抛出。

序列输入是内存中的 ASCII：`str`、`bytes`、`uint8` 数组、`types.Sequence`、它们的列表，或 `types.StringVector`（dtype 与连续性符合要求时复用输入缓冲）。结果中 `strand` 列为 `int8`：`1` 正链、`-1` 反链；MAST 合并双链的行为记为 `0`。打分预设及部分工具参数从原生层读取；Python 仍定义组合层默认值，例如 `Hits` / `Counts` 的 `pvalue=1e-4`、`strands="both"`、MAST 的 `mode="combine"` 与 MCAST 的 `statistics=True`、`seed=0`。

```python
import omicsweave as ow
from omicsweave import motif

reader = ow.open("jaspar.ffm")                       # query.MotifReader
dna = reader.alphabet
bg = motif.Background.estimate(seqs, dna, order=2)   # 二阶 Markov 背景

# 组件：两种策略合成一次扫描
fimo_like = motif.compile(reader, motif.Background.nrdb(dna), motif.Scoring.fimo())
exact = motif.compile(reader, bg, motif.Scoring.exact())
program = motif.combine(fimo_like, exact)
hits, best, counts = program.scan(
    seqs,
    motif.Hits(pvalue=1e-4, part=0),
    motif.Best(part=1),
    motif.Counts(pvalue=1e-3, part=1),
)
q = motif.stats.BenjaminiHochberg()(hits.pvalue)     # q 值是单独一步

# 打开编译结果
cm = exact.model(0, "+")
cm.null.survival, cm.score_threshold(1e-4)

# MAST 的组合 p 值换成任意 Python 函数
report = motif.mast(reader, seqs, background=motif.Background.nrdb(dna),
                    combine=lambda p: float(p.min() * len(p)))

# 工具级：输出 FIMO 风格的文本
hits = motif.fimo(reader, seqs, background=motif.Background.nrdb(dna))
print(motif.text.fimo(hits, alphabet=dna, motif_ids=reader.ids, names=names, sequences=seqs))
```

### Python 编译、扫描与结果布局

`compile(motifs, background, scoring=None, *, strands="both", alphabet=None, threads=None)` 返回 `Program`，输入参数也可写为 `motifs=...`；`combine(*programs)` 要求字母表一致，motif 序号与 part 顺序依次拼接。`Program.widths` 为 `(M,)` 的 `uint32` 数组；`part_of(i)` 返回所属 part；`model(i, strand="+")` 返回完整背景上下文的 `CompiledMotif`。`terms` 的表为 `float64`，`int_terms` 为 `int32`，每项形状为 `(n,) * len(positions)`；`null.pmf`、`null.survival` 为 `(max+1,)` 的 `float64`。模糊字母表不作为这些普通项表的额外坐标轴。

`scan(sequences, *reductions, threads=None, tile_len=None)` 至少需要一个 reduction：一个返回对应结果，多个按参数顺序返回 tuple。`Hits`、`Counts` 的 `pvalue=` 与 `score=` 互斥；均省略时为 `pvalue=1e-4`。支持 `part=` 的规约只订阅指定 part，motif 序号仍是整个程序的全局序号。`Scores` 没有 part 过滤。

下表中 `S` 为序列数、`M` 为 motif 数、`K` 为程序扫描的链数、`L_i` 为第 i 条序列长度。所有坐标是输入正向序列上的 0-based 半开区间。

| Python 结果 | 数组与 shape |
|---|---|
| `HitTable` | 每个命中一行；`sequence/motif/start/end: uint32[H]`，`strand: int8[H]`，`score/pvalue: float64[H]` |
| `BestTable` | 每个实际有分数的 `(sequence,motif)` 一行；同上但无 `end`；按序列、motif 排序，平局取更早起点、再正链 |
| `CountMatrix` | `shape=(S,M)`；`indptr: uint64[S+1]`，`motifs/counts: uint32[nnz]`；双链合计；`to_dense()` 产生新数组，`to_scipy()` 需要可选依赖 SciPy |
| `ScoreTracks` | `offsets: uint64[S+1]`、扁平 `values: float32[Σ M·K·L_i]`；若返回对象命名为 `scores`，用 `scores[i]` 取得 `(M,K,L_i)` 的视图，没有 `.tracks` 属性；无分数处 NaN；空程序通过序列长度元数据保留 `(0,0,L_i)` |
| `MastReport` | 每条序列一行，`mode="separate"` 时每链一行；`sequence: uint32[R]`，`strand: int8[R]`，`pvalue/evalue: float64[R]`，`ranked: uint64[排名行数]`；`best.indptr`、`hits.indptr` 为 `uint64[R+1]`，索引各自的扁平列 |
| `ClusterTable` | 每簇一行；`id: uint32[C]`，`sequence/start/end/hits/segment_start/segment_len: uint32[C]`，`score/excess/gc/pvalue/evalue/qvalue: float64[C]`；`threshold: float` 为全表共用的簇阈值 D，`excess=score-D`；关闭统计时 p/E/q 为 NaN，`statistics=False` |

`MastReport.best` 保存 motif、raw score、位置 p 值与 best-of 处理后的序列 p 值；`hits` 保存其非重叠排布。`ranked` 是 report 行的索引，不是序列 ID。`ClusterTable.id` 是从 1 开始的发现编号；开启统计时输出可按显著性重排，不能把行号当发现编号。直接构造 `Clusters()` 默认不计算簇统计；工具 `motif.mcast(..., statistics=True)` 默认执行统计，`statistics=False` 可关闭。

工具 `mast(motifs, sequences, *, background, letters=None, mode="combine", hit_threshold=None, evalue_threshold=None, best_of=None, combine=None, adjust=None, alphabet=None, threads=None)` 要求显式的 `background=`；`letters` 单独控制 p 值的字母频率，未提供时使用 MAST 内置频率。`mode` 支持 `"combine"`、`"separate"`、`"forward"`；字母表没有互补表时使用 `"forward"`。

工具 `mcast(motifs, sequences, *, background, pthresh=None, max_gap=None, egcost=None, threshold=None, statistics=True, seed=0, alphabet=None, threads=None)` 同样要求 `background=`，固定扫描双链；命中阈值参数名为 `pthresh`。`statistics` 可传 `True`、`False` 或自定义 `ClusterStatistics`。

结果数组独立拥有数据，下一次扫描或关闭输入读取器不会覆盖它们；`ScoreTracks[i]` 则借用当前结果的 `values`。Python 的 frozen dataclass 只约束字段重新赋值，不使 NumPy 数组自动只读。修改 `Program.model()` 导出的数组不会改动原生程序。

`compile_profiles(motifs, background, scoring=None, *, strands="both", alphabet=None)` 返回 `ProfileProgram`，默认 `ProfileScoring.glam2(alphabet)`。`scan(sequences, *, top=None, names=None, limit=None, threads=None)` 中 `top` 与 `limit` 互斥：前者保留每个 profile 的全局前 N，后者限制每条序列每链的比对数；均省略时保留全部比对。返回 `ProfileAlignments`：行列 `sequence/motif/start/end: uint32[A]`、`strand: int8[A]`、`score: float64[A]`；`column_indptr: uint64[A+1]` 索引 `matched: bool[所有列数]` 与 `offset: uint32[所有列数]`。`offset` 沿被扫描的链相对于比对起点，反链也不是正向绝对坐标；`score` 单位为 nat，`bits` 属性换算成 bit。

接口定义与列转换可定位到 [`_program.py`](../src/omicsweave/motif/_program.py)、[`_reduce.py`](../src/omicsweave/motif/_reduce.py)、[`_profile.py`](../src/omicsweave/motif/_profile.py) 与[原生绑定](../crates/python/src/motif.rs)。

## 代码结构

```text
crates/motif/src/
  background.rs   背景模型与构造、对称化、估计、MEME 背景准备
  scoring.rs      Scoring、LogScale、Quantizer、Rounding、Strands
  factor.rs       ① 估计、② 定向
  score.rs        ③ 对比（ScoreModel）、④ 量化（IntModel）
  null.rs         ⑤ 卷积、变量消元、阈值换算
  model.rs        MotifModel / StrandModel / Variant、Threshold
  lower.rs        ⑥ 分块
  scan.rs         ⑦ 链接、tile、编码、内核、并行驱动
  reduce.rs       Reducer、Demand、Match、SequenceChunk、内置规约
  program.rs      Program / Part：编译、合并与扫描入口
  stats.rs        BestOf、Combine、Adjust、TailFit 及其实现
  fimo.rs         FIMO 文本输出与 C 风格数字格式
  mast.rs         MAST：打分预设、规约、排布、相关系数、hit list
  mcast.rs        MCAST：打分预设、聚簇规约、统计、TSV 输出
  substitution.rs 替换矩阵列先验（MCAST）
  profile.rs      带 gap 的 profile：打分、局部比对、GLAM2 规约与输出
```

## 验证与性能证据

测试源码描述了待验证的断言，不代表某次构建已经运行了全部断言。当前覆盖入口包括：

| 测试 | 检查内容 |
|---|---|
| [`crates/motif/tests/integration.rs`](../crates/motif/tests/integration.rs) | 多种因子、背景阶数、逐窗口枚举、零分布、part、模糊字母以及线程/tile 交付 |
| [`crates/motif/tests/cases/`](../crates/motif/tests/cases/) | MAST 自身 context/strand 统计、线程/tile 确定性、巨大宽度/全局展开限制、prior 输出、重复编码、极端有限数值与 tiny p 值 cutoff |
| [`crates/motif/src/background.rs`](../crates/motif/src/background.rs) 的单元测试 | 频率、估计、采样上下文重启与采样边界 |
| [`src/tests/motif/`](../src/tests/motif/) | Python 读写、编译、规约、回调、统计与文本接口 |

可选黑盒用例只在对应可执行文件被显式配置时运行：

| 用例 | 环境变量 | 对照产物 |
|---|---|---|
| [`tests/external/fimo.rs`](../crates/motif/tests/external/fimo.rs) | `OMICSWEAVE_FIMO` | FIMO 文本命中 |
| [`tests/external/mast.rs`](../crates/motif/tests/external/mast.rs) | `OMICSWEAVE_MAST` | MAST hit list、排名与统计字段 |
| [`tests/external/mcast.rs`](../crates/motif/tests/external/mcast.rs) | `OMICSWEAVE_MCAST` | MCAST 聚簇文本；另设 `OMICSWEAVE_MCAST_STATS` 时包含统计对照 |
| [`tests/external/glam2.rs`](../crates/motif/tests/external/glam2.rs) | `OMICSWEAVE_GLAM2SCAN` | GLAM2Scan 比对文本 |

这些用例由 [`tests/external_tools.rs`](../crates/motif/tests/external_tools.rs) 注册，默认标记为 ignored；配置对应变量后，在有效 worker 中用 `cargo test -p omicsweave-motif --test external_tools -- --ignored` 显式运行。显式运行但缺少工具变量时会失败，不能计为对照通过。报告实际一致性时应附当前源码版本、外部工具版本、命令、启用项和日志；本页不保留无法定位原始记录的历史用时或逐字节一致性数字。

性能评估应分别量编译、已编译程序的扫描和结果转换，覆盖常见短 motif 与长跨度因子。算法降阶或局部 unchecked 不保证每种短输入都更快；比较结果需使用相同模型、序列、精度、线程与 CPU 条件，核对输出 checksum，并披露回退。内存序列的测量不代表真实 FASTA、磁盘或远程 I/O 吞吐。

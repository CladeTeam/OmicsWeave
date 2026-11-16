# Motif 文本格式

[格式总览](README.md) · [FFM](ffm.md) · [转换](../convert.md)

各数据库和工具发布 motif 的文本格式。`ow convert` 把一个文件或一组文件转成一个 [FFM](ffm.md)，之后用 Python 的 `MotifReader` 或 Rust 的 `FfmReader` / `MotifSet` 读取。所有这些文本输入也支持 gzip / BGZF 压缩；文本本身不交给 motif reader 直接查询。

```sh
ow convert JASPAR2026_CORE_vertebrates_redundant_pfms_jaspar.txt jaspar.ffm
ow convert cisbp/pwms/*.txt cisbp.ffm           # shell 展开为多个来源
ow convert hocomoco/pfm/*.pfm extra.pcm all.ffm --from pcm
```

```python
from pathlib import Path

from omicsweave import convert

convert.run(sorted(Path("cisbp/pwms").glob("*.txt")), "cisbp.ffm")
convert.run(["a.meme", "b.meme"], "all.ffm")
```

## 速览

| 格式名（`--from`） | 典型来源 | 一个文件 | 存的值 | ID / name | 转成的因子 |
|---|---|---|---|---|---|
| `meme` | MEME / STREME 输出；JASPAR、HOCOMOCO、CIS-BP 的 MEME 下载 | 多个 motif | 概率 + `nsites` | `MOTIF id alt` | 每位置一个一阶 `joint` |
| `jaspar` | JASPAR `.jaspar` | 多个 | 计数（`A [ … ]` 行） | `>ID name` | 同上，sites = 最大列和 |
| `jaspar-pfm` | JASPAR raw `.pfm`；HOCOMOCO `jaspar_format.txt` | 多个 | 计数，4 行 × W 列 | `>ID name`（无表头时用文件名） | 同上 |
| `transfac` | TRANSFAC；JASPAR、HOCOMOCO 的 TRANSFAC 导出 | 多个 | 计数或频率，列由 `P0` 声明 | `AC`（否则 `ID`）/ `ID` 行（否则 `NA`） | 同上 |
| `homer` | HOMER `.motif` | 多个 | 概率 | 表头第 2 列（ID 与 name 相同） | 同上 |
| `cisbp` | CIS-BP `pwms/*.txt` | 一个 | 概率，`Pos A C G T`（或 `U`） | 文件名 | 同上；没有行的空壳文件跳过 |
| `pcm` | HOCOMOCO `.pcm` / `.pfm`；Cluster-Buster `.cb` | 多个 | 计数或频率，W 行 × 4 列 | `>id name` | 同上 |
| `pwm` | HOCOMOCO `.pwm` | 多个 | 对数权重 | `>id` | 每位置一个一阶 `potential`，值为 `exp(权重)` |
| `dpcm` | HOCOMOCO 二核苷酸 `.dpcm` | 多个 | 计数，(W−1) 行 × 16 列 `AA AC … TT` | `>id` | 相邻位置对的二阶 `joint` |
| `dpwm` | HOCOMOCO 二核苷酸 `.dpwm` | 多个 | 对数权重 | `>id` | 相邻位置对的二阶 `potential` |
| `bamm` | BaMM `.ihbcp` | 一个 | 条件概率 | 文件名 | 每位置一个 `conditional` |
| `probound` | ProBound / MotifCentral JSON | 一个模型 | 能量系数 | 文件名 + binding mode | 单位置与位置对的 `potential` |
| `glam2` | GLAM2 `glam2.txt` | 多个比对 | 残基、缺失、插入计数 | 文件名 + 比对序号 | [带 gap 的 profile](ffm.md#带-gap-的-profile) |
| `ffm` | FFM | 多个 | 原模型数值 | ID 和 `name`，其他注释不带出 | 保留模型（用于合并、重新打包） |

当前不提供 TFFM（隐马尔可夫模型）、混合模型、IUPAC 共识序列、PBM 8-mer E-score 表的转换路由。

## 值原样保存

计数仍是计数，概率仍是概率；对数权重和 ProBound 系数存成 `exp(值)`，用于匹配时取对数。转换不归一化、不加伪计数。数值必须满足 FFM 模型约束；例如 `exp(值)` 溢出为无穷时不能写入。

计数与概率保留的信息不同，压缩后的大小也依赖数值精度和分布；转换不会从概率自动反推计数。背景、阈值及匹配时的伪计数不写进模型。

未显式提供 sites 的 `pcm`、HOMER、CIS-BP、TRANSFAC，以及 `dpcm`，共用计数推断：取各位置表的数值和，最大和大于 `1 + 1e-3` 时将其设为 sites，否则不设。JASPAR / raw PFM 则显式取最大列和；MEME 优先保留正的 `nsites`。

## 注释

每个 motif 只保留一个注释键 `name`：MEME 的 alt ID、JASPAR / HOCOMOCO 表头中 ID 之后的部分、TRANSFAC 以 `AC` 为 ID 时的 `ID` 行（否则 `NA`）、HOMER 的名字列；没有时为空。URL、E 值、HOMER 阈值、物种等字段不保留，它们要么是扫描参数，要么属于数据库的元数据表。

下文用文件名生成 ID 时，取去掉目录、可选 `.gz` 后缀及最后一个扩展名后的部分。例如 `pwms/TF1.txt.gz` 的 ID 为 `TF1`。

## 多个文件

来源可以是一个文件，或一组文件：

- 来源序列不能为空，按给出的顺序合并；相同路径给出多次只处理一次；
- 所有文件必须是同一种格式，且字母表定义相同，包括名称、符号顺序、互补和别名（DNA 与 RNA 的 CIS-BP 文件需分开转换）；
- motif ID 在整个集合中必须唯一，重复即报错；
- 任一文件不存在时报错。

OmicsWeave 不展开通配符：CLI 中由 shell 把 `pwms/*.txt` 展开成多个参数（最后一个参数是目标文件）；Python 中传路径列表。

多源转换的 `plan.source` / `report.source` 只表示第一个来源。空的 glob 在 Python 中得到空列表时会报错。显式格式对象或 `--from` 会指定所有来源的解析器，之后仍会检查每个来源的内容。

## 内存与缓存

当前转换按文件整体解析，文本和该文件的模型会占用内存；FFM 来源按 block 解码一次，装入 resident 集合，再通过借用视图写入。它不依赖 reader 的隐式跨调用缓存，也不为转换创建共享 `BlockCache`。压缩后的输出块保留到最终写出，整个路由不 spill 到磁盘，`--tmp-dir` 不能把它变成外存转换。

计划中的固定内存为源文件总字节数的四倍估算。压缩输入的展开体积、模型元数据和写入期间的缓冲会使实际内存超过这个估算。运行时限制展开文本长度并给 FFM 解码设置分配上限，但 `memory` 不是整次转换的 RSS 硬限制；大集合需预留实际内存，详见 [资源预算](../convert.md#资源预算试运行和输出)。

## 格式识别

不指定 `--from` 时按内容识别，每个来源都识别一次，结果必须一致。文本探测最多读取开头 64 KiB（压缩文本先展开），不以来源扩展名决定格式：

| 特征 | 格式 |
|---|---|
| 首行 `MEME version` | `meme` |
| 首行 `Pos A C G T` | `cisbp` |
| 两字母行代码（`AC`、`ID`、`XX` …）后跟 `P0` / `PO` 表头 | `transfac` |
| JSON 对象且含 `"coefficients"` | `probound` |
| 空行分组的数字行，各组的列数依次为 `A, A², A³, …`，且组间阶数相同 | `bamm` |
| 首行 `GLAM2: Gapped Local Alignment of Motifs` | `glam2` |
| `>` 表头后的某行同时含 `[` 和 `]` | `jaspar`；探测阶段不要求行首为 `A`，实际解析再校验字母和数值 |
| `>` 表头第 3 列是数字、每行 4 个数 | `homer` |
| `>` 表头下每行 16 个数 | `dpcm`，有负数时 `dpwm` |
| `>` 表头下每个 motif 都是 4 行、各行等宽且无负数，并排除上面的布局 | `jaspar-pfm` |
| `>` 表头下每行 4 个数 | `pcm`，有负数时 `pwm` |
| `>` 表头下是序列字母 | FASTA |

4 × 4 的方阵既可能是 4 个位置的 `pcm`，也可能是宽 4 的 `jaspar-pfm`，无法区分，此时报错并要求用 `--from` 指定。所有非负值的权重矩阵也可能被识别为计数格式；已知来源是权重时应显式使用 `pwm` / `dpwm`。BaMM 可按上述数字结构自动识别；最后一个数字块可以直接在文件末尾结束，不要求尾随空行。探测非完整前缀时，若已有多个块且没有以双换行结束，会先丢弃可能截断的最后一块。探测前缀不足或布局有歧义时用 `--from bamm`。格式探测不是完整输入校验。

## 各格式细节

### MEME

读取 `ALPHABET=` 短格式（`ACGT` 为 DNA、`ACGU` 为 RNA、20 种氨基酸为蛋白）或 `ALPHABET "name" … END ALPHABET` 自定义字母表块：核心符号成为 FFM 字母表的符号，互补对成为互补表，只对应一个核心符号的别名成为别名，模糊符号（如 `N = ACGT`）映射为 invalid。没有字母表声明时，下一张 `letter-probability matrix` 表头的 `alength=20` 或背景列出 20 种符号才会猜为蛋白，不直接数矩阵行的列数；背景有 `U` 而无 `T` 时猜为 RNA，否则为 DNA。

只读 `letter-probability matrix`：`alength=` 必须等于字母表大小；`w=` 必须是非负 32 位整数，提供时行数必须相符，否则读到第一行非数字为止；正的 `nsites=` 成为 sites。没有矩阵或矩阵为空的 motif 跳过并警告。背景可用于猜测字母表，但与 `strands:`、`log-odds matrix`、`URL` 一样不作为模型字段保留。

### JASPAR 与 raw PFM

`.jaspar` 的行字母可以任意顺序；raw PFM 按 A、C、G、T 顺序。两者都是计数，sites 取各列之和的最大值。

### TRANSFAC

列由 `P0`（或 `PO`）行声明，可以任意顺序，含 `U` 时为 RNA；位置行末尾的共识字母忽略。sites 使用上面的通用计数推断规则。

### HOCOMOCO 与 Cluster-Buster

`.pcm`、`.pfm`、Cluster-Buster 同为「每位置一行」，由 `pcm` 读取；以 `#` 开头的行忽略。逐 motif 的文件一次转换，例如 `ow convert H14CORE/pfm/*.pfm h14.ffm`。

`pwm` 与 `dpwm` 按自然对数权重解释，写入 `exp(权重)`，转换不再添加背景或伪计数。二核苷酸输入每个 motif 有 (W−1) 行，每行 16 列，依次为 `AA AC AG AT CA … TT`（前一个位置变化慢）。名为 `jaspar_format.txt` 的下载也需核对布局：4 行无标签计数由 `jaspar-pfm` 读取，`A [ … ]` 形式则由 `jaspar` 读取。

### HOMER

表头以 tab 分隔：共识序列、名字、对数似然阈值，以及可选的统计列；只用名字列。数值行可以是概率，也可以是未归一化的值。

### CIS-BP

每个文件一个 motif，表头为 `Pos` 加 A、C、G、T（或 U），字母列可以重排。首列位置编号不参与重新排序，矩阵按输入行顺序保存。只有表头、没有矩阵行的文件跳过并警告；大量跳过时合并为汇总警告。

### BaMM

BaMMmotif2 的 `.ihbcp`：没有表头，每个位置一组 `K + 1` 行，组间用空行分隔，最后一组可直接在 EOF 结束；第 `k` 行是 `A^(k+1)` 个 `k` 阶条件概率，当前碱基变化最快、最早的上下文碱基最慢，与 FFM 因子表的排列相同。BaMM 在 motif 开头把低阶概率平铺成整行；FFM 保留真实阶数：位置 `j` 的因子阶数为 `min(j, K)`，覆盖位置 `j − k … j`。字母表由第一行的宽度判断：4 列为 `ACGT`，6 列为 BaMM 的 `ACGTMH`；5 列可能是 `ACGTM` 或 `ACGTH`，文件中无法区分，因此拒绝转换。解析器不按扩展名或概率含义区分联合概率 `.ihbp` 与条件概率 `.ihbcp`：布局相同的 `.ihbp` 也可能识别为 `bamm` 并按条件概率保存，不会自动换算。输入必须确实符合上述条件模型语义；背景文件 `.hbcp` / `.hbp` 没有单独的转换路由。

### ProBound

MotifCentral 的 fit 是一个 JSON 对象；ProBound 的输出是 JSON lines，取最后一行（最终模型）。字母表来自 `modelSettings.letterOrder`（`ACGT` 为 DNA，其他按其顺序建立，并按 `letterComplement` 设互补）。每个 `size` 大于 0 的 binding mode 是一个 motif：

- `mononucleotide[n·x + letter]` 成为位置 `x` 的一阶 `potential`；
- `dinucleotide[d − 1][n²·x + n·a + b]` 成为位置对 `(x, x + d)` 的二阶 `potential`；

值为 `exp(系数)`，与 ProBound 的相对亲和力 `exp(Σ 系数)` 一致。`activity`、`positionBias`、`modifications` 和 mode 之间的相互作用不属于单个 mode 的模型，不转换。ID 为文件名；有多个特异 mode 时为 `文件名.i`，`i` 是原 `bindingModes` 数组中从 0 起的索引，跳过非特异 mode 后不重新编号。

### GLAM2

读取 GLAM2 的文本输出 `glam2.txt`（`glam2.html` 与 `glam2.meme` 不是转换源；后者只有概率，丢失了 gap）。文件中每个比对以 `Score: … Columns: W  Sequences: N` 开头，以列表结束：列行是各残基的计数、`Del`（该列被缺失的序列数）和列得分，相邻列之间的行是 `Ins`（插入的残基总数）和得分。每个比对成为一个[带 gap 的 profile](ffm.md#带-gap-的-profile)：

- 每列一个一阶 `joint` 因子，表为该列的残基计数；
- `deletions` 为 `W` 个 `Del`，`insertions` 为 `W − 1` 个 `Ins`；
- sites 为 `Sequences:`（参与比对的序列数）。

得分由计数推出，不保存。字母表来自列表表头：`a c g t` 为 DNA，20 种氨基酸为蛋白，其他按 GLAM2 所用字母表文件的符号建立（同义字母和通配符不在输出中，无法还原）。某列在所有序列中都被缺失时，该列的残基计数全为 0，原样保存。没有列的比对跳过并警告。ID 为文件名；文件中有多个比对时为 `文件名.k`（`k` 从 1 开始，与 GLAM2 的编号一致）。GLAM2 的输出都叫 `glam2.txt`，合并多次运行的结果时需先改名，否则 ID 重复。

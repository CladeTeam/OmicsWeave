# MAF

[格式总览](README.md) · [转换](../convert.md) · [查询](../query.md) · [bigMaf](bigmaf.md)

MAF（Multiple Alignment Format）是保存**多物种基因组比对**的文本格式：每个比对块把多个物种中同源的一段序列逐列对齐，块内每个物种一行。UCSC 的 multiz 比对（如 hg38 的 `multiz100way`、`multiz470way`）、Ensembl 的 EPO 比对、Cactus 输出的 `.maf` 都使用它。phyloP、phastCons 等保守性分数就是从这些比对计算出来的；直接使用比对本身作为输入的模型（如 GPN-MSA）则需要按区域取出各物种的比对列。

在 OmicsWeave 中，MAF 是**转换源**：`ow convert` 把它转成 [bigMaf](bigmaf.md)（一种 bigBed），再由 `MsaReader` 按区域读取为 `[物种, 长度]` 矩阵。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 多物种比对块：每个物种一行带 gap 的序列 |
| 常见后缀 | `.maf`、`.maf.gz` |
| 存储 | 空格分隔文本，可用 gzip / bgzip 压缩 |
| 坐标 | 0-based；负链行的起点从反向互补序列的 5' 端起算 |
| 随机访问 | 不支持；转为 bigMaf 后支持 |
| 典型来源 | UCSC `goldenPath/<genome>/multiz*way/maf/`、Ensembl、Cactus（`cactus-hal2maf`） |
| OmicsWeave | 格式识别、转换为 bigMaf |

## 文件结构

```text
##maf version=1 scoring=roast.v3.3
a score=25976.000000
s hg38.chr22             10510000 2 +  50818468 GA
s panTro4.chrUn_GL393523    88624 2 +    405060 GA
i panTro4.chrUn_GL393523 N 0 C 0
s ponAbe2.chrUn          66603503 2 -  72422247 GA
e mm10.chr16             25201113 87 + 98207768 I

a score=3021.000000
s hg38.chr22 10510002 5 + 50818468 TC-CAG
s rheMac3.chr7 24246326 6 + 170124641 TCACAG
```

文件以 `##maf version=1` 开头，之后是由空行分隔的段落；`#` 开头的行是注释。每个比对块以 `a` 行开始：

| 行 | 内容 |
|---|---|
| `a score=... ` | 块的开始，`key=value` 形式的属性，通常有 `score`、`pass` |
| `s src start size strand srcSize text` | 一个物种的比对行 |
| `i src leftStatus leftCount rightStatus rightCount` | 前一个 component 与相邻块的衔接状态（通常跟在 `s` 后；当前解析器也接受跟在 `e` 后） |
| `e src start size strand srcSize status` | 该物种在此处没有可比对的序列，但两侧有比对（`I` 插入、`C` 连续、`M` 缺失数据等） |
| `q src quality` | 前一个 component 的质量文本（通常跟在 `s` 后；当前解析器也接受跟在 `e` 后） |
| `r type size id` | 区域定义，少见 |

`s` 行各字段：

| 字段 | 含义 |
|---|---|
| `src` | `物种.染色体`，如 `hg38.chr22`；物种名与染色体名之间用第一个点分隔 |
| `start` | 0-based 起点。负链时相对反向互补序列计算 |
| `size` | 比对中该物种的碱基数（不含 `-`） |
| `strand` | `+` 或 `-` |
| `srcSize` | 整条染色体的长度 |
| `text` | 比对文本，`-` 为 gap；一个块内所有行的列数相同 |

## 关键概念

### 参考物种

multiz 风格的 MAF 以一个物种为**参考**：每个块的第一行都是参考物种，且在正链上，块按参考坐标排列。其他物种可以出现在任意链、任意染色体上。负链行的 forward 区间为 `[srcSize − start − size, srcSize − start)`；其 `text` 已按当前比对块方向书写，投影时不再因该行的 `strand` 单独互补。bigMaf 用参考行的坐标作为每个块的区间，查询也按参考坐标进行。

### 参考上的 gap

块内某列参考为 `-`，表示其他物种在此有参考没有的插入。投影到参考坐标时，这些列没有对应的参考位置，会被丢弃；其他物种行中的 `-` 则表示该物种相对参考的缺失。

### 没有比对的物种

一个块不一定包含所有物种：某物种在这个区域没有可比对的序列时，可能有一行 `e`，也可能完全没有。没有任何块覆盖的参考位置也存在。读取时这两种情况都用填充字符表示。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 文件对象 | `Maf`（纯文本）、`GzipMaf`（gzip 或 bgzip） |
| 格式识别 | 内容以 `##maf` 开头 |
| 区域查询 | 不直接支持；`MsaReader`、`query.open` 会报错并给出转换命令 |
| 转为 bigMaf | 每个 `a` 块一条 bed3+1 记录，采用 `mafToBigMaf` 的文本布局，见 [bigMaf](bigmaf.md) |

```sh
ow convert chr22.maf.gz chr22.bigmaf
ow convert multiz.maf.gz multiz.bigmaf --ref-species hg38 -t 16 -m 8G
```

```python
from omicsweave import convert

convert.run("chr22.maf.gz", "chr22.bigmaf", ref_species="hg38", threads=16)
```

当前转换的检查和限制：

- `##maf` 头必须存在且 `version=` 解析为非零整数；`version=0` 与缺少有效版本一样报 `no version on the ##maf line`。注释、`##maf` 头和不以 `a` 开头的段落被跳过。
- `s` / `e` 行必须是 7 个字段；整数按有符号 32 位范围解析，`start`、`srcSize` 不能为负。`s` 行还检查 `size >= 0` 且 `start + size <= srcSize`，同一块的 `s` 行文本列数必须相同。
- `i`、`q` 行的 `src` 必须与它们之前的 component 相同，该 component 可以是 `s` 或 `e`；输出 `e` 时不会另写其 `q` / `i` 行。
- 参考物种默认取第一个块首个 component 的 `src` 在第一个点之前的部分，也可用 `ref_species` 指定。它必须是每个块的第一个 component，且在 `+` 链上；参考行的 `srcSize` 给出染色体长度，同一染色体的声明必须一致。
- `chrom_sizes` 可选，提供时必须与参考行的 `srcSize` 一致，不会改写输入坐标或截断块。
- 输入不要求按参考坐标排序；块在参考上重叠会报错。

语法错误带有行号。转换保留 component 的比对文本，不过滤物种、不投影；`a` 行只保留 `score`，其他属性丢弃，输出细节见 [bigMaf](bigmaf.md)。它不替代完整的 MAF 校验：例如转换阶段不逐行核对 `size` 是否等于文本的非 gap 数量；`MsaReader` 查询时另行检查所有读到的 `s` 行：七个字段、UTF-8 的 `src`、可打印 ASCII 的比对文本、非负整数坐标、合法链、`start + size ≤ srcSize`、各行列数一致及非 gap 数量等于 `size`；即使该物种未被选择也会检查。首个 `s` 行还必须在正链上，且坐标、染色体后缀与 `srcSize` 分别匹配 bigBed 区间和染色体表。转换成功不代表所有这些查询约束已经检查；查询不解析 `a` 的得分或 `e/i/q` 的附加语义。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| `not with the reference ...` | 有的块不以参考物种开头；确认 `--ref-species`，或先用 `mafSpeciesSubset`、`maf-swap` 等工具调整 |
| `the reference is not on the + strand` | 参考行在负链；bigMaf 要求参考在正链 |
| `records overlap` | 块在参考上重叠，常见于未经 multiz 合并的两两比对 |
| `N long here but M in an earlier block` 或 `chromosome ... is N long here but M ...` | 同一参考染色体的 `srcSize` 前后不一致；同一解析批次与跨批次汇合的报错文案不同 |
| 查询时报行的 bases / columns 不一致 | 检查每个 `s` 行的 `size` 是否等于非 `-` 字节数，以及所有行的比对文本是否同长；未选中的物种也不能携带畸形行 |
| 物种名含点 | 转换时可显式指定含一个点的参考名，如 `ref_species="hg19.1"`；`MsaReader` 当前仍按第一个点分组，行名会是 `hg19` |

## 相关格式

- [bigMaf](bigmaf.md)：MAF 转换后的可查询格式。
- [bigBed](bigbed.md)：bigMaf 所用的容器。
- [bigWig](bigwig.md)：phyloP、phastCons 等保守性分数的常见格式，可直接用 `TrackReader` 读取。

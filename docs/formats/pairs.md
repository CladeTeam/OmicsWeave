# 4DN pairs

[格式总览](README.md) · [转换](../convert.md) · [定量](../quant.md) · [.hic](hic.md) · [Juicer short](juicer-short.md) · [cool](cool.md)

4DN pairs（`.pairs`）是 4D Nucleome 项目制定的**接触对文本格式**：Hi-C、Micro-C 等实验比对后，每一对读段的两端位置写成一行。它位于 BAM 与接触矩阵之间：从 BAM 解析出 pairs，再按某个分辨率分箱成 [`.hic`](hic.md) 或 [cool](cool.md) 矩阵。

OmicsWeave 把 pairs 作为 `.hic` 的转换源，可一次生成多种分辨率；`quant` 从 Hi-C BAM 计算接触时也会输出 pairs 文件。pairs 本身不做区域查询。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 每行一个接触：两端的染色体、位置、链及其他字段 |
| 常见后缀 | `.pairs`、`.pairs.gz`（通常 bgzip） |
| 存储 | tab 分隔文本，带 `#` header |
| 坐标 | `pos1` / `pos2` 为 1-based 位置 |
| 随机访问 | 可借助 pairix 索引（`.px2`）；OmicsWeave 不需要 |
| 典型来源 | pairtools、distiller、Juicer 工具链、`ow quant --hic` |
| OmicsWeave | 识别、转换为 `.hic`；quant 输出 |

## 文件结构

```text
## pairs format v1.0.0
#sorted: chr1-chr2-pos1-pos2
#shape: upper triangle
#genome_assembly: hg38
#chromsize: chr1 248956422
#chromsize: chr2 242193529
#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 pair_type
read1	chr1	10001	chr1	50001	+	-	UU
read2	chr1	20500	chr2	1200	-	+	UU
```

### header

header 由 `#` 开头的行组成，第一行必须是 `## pairs format v1.0...`：

| header 行 | 含义 |
|---|---|
| `## pairs format v1.0` | 格式声明 |
| `#sorted:` | 排序方式，例如 `chr1-chr2-pos1-pos2` |
| `#shape:` | `upper triangle` 表示两端已按染色体和位置排序，使 end1 不在 end2 之后 |
| `#genome_assembly:` | genome 名称 |
| `#chromsize:` | 一行一条染色体：名称与长度，顺序即染色体顺序 |
| `#columns:` | 数据列名 |
| `#samheader:` | 可选，保存来源 BAM 的 header |

### 数据列

前七列是标准列：

| 列 | 含义 |
|---|---|
| `readID` | 读段名；可为 `.` |
| `chr1`、`pos1` | 第一端的染色体和 1-based 位置 |
| `chr2`、`pos2` | 第二端的染色体和位置 |
| `strand1`、`strand2` | 两端的方向，`+` / `-` |

之后可以有任意附加列，常见的有：

- `pair_type`：pairtools 的分类，如 `UU`（两端唯一比对）、`UR` / `RU`（一端被“救回”）、`NN`、`MM`、`WW`（walks）等。
- `mapq1`、`mapq2`：两端的比对质量。
- `frag1`、`frag2`：所属酶切片段。

位置通常是读段 5′ 端（连接点附近）在参考上的坐标，具体取法由生成工具决定。

### pairix 索引

bgzip 压缩并排序的 pairs 可以用 pairix 建 `.px2` 索引，按两个区域查询接触对。OmicsWeave 转换时顺序读取全文件，不读取也不需要 pairix 索引。

## 关键概念

### 接触对与接触矩阵

pairs 中每行是一次观测；矩阵中每个格子是落在两个 bin 中的观测数。分箱规则是 `bin = position // resolution`，因此 pairs 可以生成任意分辨率，而矩阵一旦分箱就不能变得更细。保留 pairs 意味着以后可以重新选择分辨率或过滤条件。

### 过滤在哪一步完成

PCR 重复、低 MAPQ、`pair_type` 不合格、同一酶切片段内的自连接等，通常在生成 pairs 时由 pairtools 等工具过滤或标注。分箱步骤只计数，不再判断这些字段。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 格式识别 | `formats.detect` 按 `## pairs format` 识别为 `Pairs` / `GzipPairs` |
| 压缩 | 纯文本、gzip、BGZF |
| 转换 | pairs → `.hic`（v9），一次生成多种 BP 分辨率 |
| 染色体 | 默认取 `#chromsize:`；可用 `chrom_sizes` 覆盖 |
| 列 | 按 `#columns:` 定位 `chr1`、`pos1`、`chr2`、`pos2`（也接受 `position1` / `position2`）；无该行时按标准七列 |
| genome 名称 | 取 `#genome_assembly:`，缺失时写 `unknown` |
| quant 输出 | `ow quant --hic` 写出 `hic/contacts.pairs.gz` |

转换的计数规则：

- 每行计一次接触，不读取 MAPQ、`pair_type`、重复或片段字段，不做任何过滤。
- 位置按文本值直接分箱：`bin = position // resolution`，不先减一；超出染色体长度的位置先截断为染色体长度再分箱（与 Juicer `pre` 相同）。
- 两端先按染色体顺序排列，同一染色体内的接触放入上三角，所以 `(a, b)` 和 `(b, a)` 计入同一个格子。
- 染色体表之外的行被跳过，并按染色体汇总写入转换 warning。

当前只输出 BP 分辨率，不生成 FRAG 分辨率、expected 或归一化向量；不直接对 pairs 做区域查询。

## 转换

```sh
ow convert contacts.pairs.gz contacts.hic --resolutions 5000,25000 -t 4 -m 2G
```

```python
import warnings
from omicsweave import convert
from omicsweave.query import ContactReader

report = convert.run("contacts.pairs.gz", "contacts.hic", resolutions=[5000, 25000])
print(report.records)   # 计入矩阵的接触行数（不含被跳过的行）
print(report.warnings)  # 例如被跳过的 chrUn_* 行数

with ContactReader("contacts.hic") as contacts:
    print(contacts.resolutions())  # (25000, 5000)
```

- 未指定 `resolutions` 时使用 Juicer `pre` 的九档默认分辨率：2.5 Mb、1 Mb、500 kb、250 kb、100 kb、50 kb、25 kb、10 kb、5 kb。
- 较粗的分辨率由最细分辨率的像素相加得到，因此每一档都必须是最细分辨率的整数倍，例如 `5000,25000,100000` 可以，`5000,12000` 不行。
- header 没有 `#chromsize:` 时必须提供 `chrom_sizes` / `--chrom-sizes`（两列长度表、`.fai` 或 `.2bit`）。提供 `chrom_sizes` 时，它决定染色体集合和顺序，不在其中的行被跳过。
- 输入行不需要排序。转换在缓冲预算内排序，超出时写临时文件；临时文件位置用 `tmp_dir` / `--tmp-dir` 指定。

只想保留主要染色体时，可以给一个只含 `chr1`–`chrX` 的长度表，其余染色体的行会被跳过并计入 warning。

### 需要过滤时

转换不看 `pair_type` 或 MAPQ。需要只保留 `UU` 或高 MAPQ 的接触时，先在 pairs 上过滤：

```sh
pairtools select '(pair_type == "UU") and (mapq1 >= 30) and (mapq2 >= 30)' \
  contacts.pairs.gz -o filtered.pairs.gz
ow convert filtered.pairs.gz filtered.hic --resolutions 5000,25000
```

## 从 BAM 生成 pairs

```sh
ow quant hic.bam -o hic-results --hic --hic-resolutions 5000,25000
ow quant hic.bam -o hic-results --hic --hic-resolutions none   # 只输出 pairs
```

quant 取两个 mate 的 5′ 位置，做 MAPQ 过滤和精确去重，写出 `pair_type` 为 `UU` 的 pairs，不生成 pairix 索引。它只使用主比对，不解析 walks 或嵌合多段比对，与 pairtools 的完整流程不同。规则见 [quant](../quant.md)。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 报没有 `#chromsize` | 用 `--chrom-sizes` 提供染色体长度 |
| warning 说跳过了一些行 | 这些行的染色体不在染色体表中，通常是 decoy、alt 或 unplaced contig |
| 报分辨率不是最细分辨率的倍数 | 调整分辨率组合，使每一档都能被最细分辨率整除 |
| 结果中仍有低质量接触 | 转换不过滤；先用 pairtools 等工具筛选 pairs |
| 报 `#columns` 缺少某列 | 列名需要包含 `chr1`、`pos1`、`chr2`、`pos2` |

## 相关格式

- [Juicer short](juicer-short.md)：Juicer 的另一种接触对文本，没有 header。
- [.hic](hic.md)：转换目标，按窗口查询接触矩阵。
- [cool](cool.md)、[mcool](mcool.md)：cooler 的接触矩阵，可由 pairs 分箱得到。
- [BAM](bam.md)：quant 从 Hi-C BAM 生成 pairs。

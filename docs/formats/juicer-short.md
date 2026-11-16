# Juicer short

[格式总览](README.md) · [转换](../convert.md) · [.hic](hic.md) · [4DN pairs](pairs.md)

Juicer short 是 Juicer 工具链 `pre` 命令接受的**接触对文本格式**之一：每行一个接触，只有八列（或带 score 的九列），没有 header。它是从 Juicer 流程的 `merged_nodups.txt` 精简出来的最小输入，也常被其他工具当作生成 `.hic` 的中间格式。

OmicsWeave 把 Juicer short 作为 [`.hic`](hic.md) 的转换源。由于文件没有染色体长度信息，转换时必须提供 `chrom_sizes`。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 每行一个接触，可选带权重 |
| 常见后缀 | 无固定后缀，常见 `.txt`、`.short`、`.txt.gz` |
| 存储 | 空白分隔文本，没有 header |
| 坐标 | 位置按文本值使用；Juicer 流程中通常为 1-based |
| 随机访问 | 不支持 |
| 典型来源 | Juicer 流程、从 pairs 或 BAM 自行导出 |
| OmicsWeave | 识别、转换为 `.hic`（需要 `chrom_sizes`） |

## 文件结构

```text
str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 [score]
```

| 列 | 含义 |
|---|---|
| `str1` / `str2` | 两端读段方向：`0` 表示正向，非 0（Juicer 流程中常为 `16`）表示反向 |
| `chr1` / `chr2` | 两端染色体 |
| `pos1` / `pos2` | 两端位置 |
| `frag1` / `frag2` | 两端所属酶切片段编号；不关心片段时可写 `0`，或令两端不同 |
| `score` | 可选第九列：这条接触的权重，可以是小数 |

例如，一条两端都在 chr1 的接触，以及一条带权重 2.5 的接触：

```text
0 chr1 10001 0 16 chr1 50001 0
0 chr1 10001 0 16 chr2 20001 1 2.5
```

Juicer 的 `pre` 要求同一文件中的行要么都是八列，要么都是九列（short with score）。没有 score 时每行计一次；有 score 时这一行贡献 score 个接触。

Juicer 还有 medium 和 long 两种更宽的格式（包含读段名、MAPQ、CIGAR、序列等）。OmicsWeave 只接受 short 及其带 score 的形式。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 格式识别 | `formats.detect`：所有数据行都符合八 / 九列 short 结构时识别为 `JuicerShort` / `GzipJuicerShort` |
| 压缩 | 纯文本、gzip、BGZF |
| 分隔符 | 空格或 tab |
| 转换 | Juicer short → `.hic`（v9），一次生成多种 BP 分辨率 |
| 染色体 | 必须由 `chrom_sizes` 提供（长度表、`.fai` 或 `.2bit`），同时决定染色体顺序 |
| score | 第一条数据行为九列时，全部行按 score 计数，支持小数 |

转换规则：

- 每行计一次接触（有 score 时计 score），只使用染色体、位置和 score；不根据方向或 fragment 字段过滤。
- 位置直接分箱：`bin = position // resolution`，不先减一；超出染色体长度的位置先截断为染色体长度再分箱（与 Juicer `pre` 相同）。
- 同一染色体内的接触放入上三角，`(a, b)` 与 `(b, a)` 计入同一格子。
- 染色体表之外的行被跳过，并写入转换 warning。
- 带小数 score 时，各分辨率的求和顺序固定，结果与线程数和内存预算无关。

只输出 BP 分辨率，不生成 FRAG 分辨率、expected 或归一化向量。

## 转换

```sh
ow convert contacts.short contacts.hic --from juicer-short \
  --chrom-sizes genome.chrom.sizes --resolutions 5000,25000
```

```python
from omicsweave import convert
from omicsweave.formats import JuicerShort

convert.run(
    JuicerShort("contacts.txt"), "contacts.hic",
    chrom_sizes="genome.chrom.sizes", resolutions=[5000, 25000],
)
```

内容符合 short 结构时会被自动识别；文件名没有特征后缀、或者希望明确表达意图时，用 `JuicerShort(...)` / `GzipJuicerShort(...)` 或 CLI 的 `--from juicer-short` 显式指定。

`resolutions` 默认是 Juicer `pre` 的九档 2.5 Mb 到 5 kb；较粗的每一档必须是最细分辨率的整数倍。转换结果按窗口查询，见 [.hic](hic.md)。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 被识别为 BED 或其他格式 | 某些行不符合八 / 九列结构；检查文件，或显式指定 `--from juicer-short` |
| 报缺少染色体长度 | short 没有 header，必须提供 `--chrom-sizes` |
| 某些行被跳过 | 它们的染色体不在 `chrom_sizes` 中，数量见 warning |
| 有的行八列、有的行九列 | 同一文件中列数必须一致 |
| 需要按 MAPQ 或片段过滤 | short 不含 MAPQ；在生成 short 之前过滤 |

## 相关格式

- [4DN pairs](pairs.md)：带 header 的接触对文本，可直接读取染色体长度。
- [.hic](hic.md)：转换目标。

# HiC-Pro valid pairs

[格式总览](README.md) · [转换](../convert.md) · [.hic](hic.md) · [4DN pairs](pairs.md)

HiC-Pro 的 `.validPairs` / `.allValidPairs` 保存每行一个接触对，没有表头。OmicsWeave 接受纯文本、gzip 和 BGZF，将它们转换为 `.hic` v9；必须提供 `chrom_sizes`，同时指定染色体名称、长度和顺序。

## 文件结构与坐标

```text
readID chr1 pos1 strand1 chr2 pos2 strand2 fragment_size res_frag1 res_frag2 mapq1 mapq2 [allele]
```

字段以 tab 分隔，共 12 列，等位特异分析可多一列 `allele`；读取也接受空白分隔，允许行尾多余的 tab。第 4、7 列方向必须是 `+` / `-`，第 3、6 列是 `0..=2**64-1` 内的非负整数位置；超出这个整数范围会报错，超出染色体长度则按下面的规则截断。内容识别按这些结构条件区分数字方向的 Juicer short 和带 `## pairs format` 表头的 4DN pairs；它不检查 fragment / MAPQ / allele 的语义，也不依据来源扩展名决定格式。

HiC-Pro 的位置是 **1-based 读段 5′ 端**。其 [`mapped_2hic_fragments.py`](https://github.com/nservant/HiC-Pro/blob/master/scripts/mapped_2hic_fragments.py) 中，正链位置来自 `reference_start + 1`，负链来自 `reference_start + alen`；负链使用 `alen` 计算位置，OmicsWeave 不重新推断 CIGAR 或读段端点。[`hicpro2juicebox.sh`](https://github.com/nservant/HiC-Pro/blob/master/bin/utils/hicpro2juicebox.sh) 原样把两端位置交给 Juicer `pre`。

转换遵循同一规则：**直接按输入值分箱**，`bin = position // resolution`，不减一。例如位置 `5000` 在 5 kb 分辨率落入起点为 `5000` 的 bin。超出染色体长度的位置截断为该长度。查询 `.hic` 的窗口仍然使用 OmicsWeave 的 0-based 半开区间。

## 转换

```sh
ow convert sample.allValidPairs sample.hic --chrom-sizes genome.sizes
ow convert sample.validPairs.gz sample.hic --chrom-sizes genome.sizes \
  --resolutions 5000,25000
```

```python
from omicsweave import convert
from omicsweave.formats import HicProPairs

convert.run(
    HicProPairs("sample.allValidPairs"), "sample.hic",
    chrom_sizes="genome.sizes", resolutions=[5000, 25000],
)
```

内容符合结构时自动识别；也可用 `HicProPairs` / `GzipHicProPairs` 或 `--from hicpro-pairs` / `--from gzip-hicpro-pairs` 指定。

每行计数一次，不按 MAPQ、酶切片段、重复读段或 allele 过滤。同一染色体内的接触合并到上三角；跨染色体接触按染色体表顺序归并。表中没有的染色体所在行被跳过，warning 按染色体汇总数量。不做染色体别名映射。

输入不需要预排序；转换复用 pairs 的有界排序与磁盘 spill 流程，输出不随线程数或 spill 划分变化。`report.records` 是计入矩阵的行数，不包含未知染色体的行，也不是聚合后的非零像素数。HiC-Pro 没有 genome assembly header，输出的 genome 标识为 `unknown`。

默认输出 Juicer `pre` 的九档 BP 分辨率（2.5 Mb 到 5 kb）。较粗分辨率必须是最细分辨率的整数倍。只写原始计数，不生成 FRAG 分辨率、expected、归一化向量或 Juicer 的 `All` 染色体。与其他工具比较时应选择 `NONE` 归一化，并核对对方的 MAPQ、去重和同片段过滤选项。

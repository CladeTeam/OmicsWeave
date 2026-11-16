# mcool

[格式总览](README.md) · [转换](../convert.md) · [cool](cool.md) · [.hic](hic.md)

`.mcool` 是**多分辨率的 cooler 文件**：在一个 HDF5 文件中，为同一份接触数据保存多种 bin 宽度的矩阵，例如 1 kb、5 kb、10 kb、…、1 Mb。HiGlass 等浏览器用它在缩放时切换分辨率，4DN 数据门户也以 mcool 发布大量 Hi-C 数据。

OmicsWeave 把 mcool 转换为同样是多分辨率的 [`.hic`](hic.md)，可以转换全部分辨率，也可以只挑选其中几种。转换只搬运原始计数，不计算新的分辨率。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 同一接触矩阵的多种分辨率 |
| 常见后缀 | `.mcool` |
| 存储 | HDF5；每个分辨率是一个完整的 cooler 组 |
| 坐标 | 与 [cool](cool.md) 相同，0-based 半开碱基坐标 |
| 随机访问 | 格式本身支持；OmicsWeave 用于转换，不直接查询 |
| 典型来源 | `cooler zoomify`、distiller、4DN 数据门户 |
| OmicsWeave | 识别、转换为 `.hic`（全部或所选分辨率） |

## 文件结构

```text
/
└── resolutions/
    ├── 1000/      一个完整的 cooler：chroms/ bins/ pixels/ indexes/
    ├── 5000/
    ├── 10000/
    └── ...
```

每个 `/resolutions/<bin_size>` 组的内部结构与单个 [`.cool`](cool.md) 文件的根组完全相同，有自己的 bins、pixels、索引和属性。各分辨率描述的是同一个基因组，染色体名和长度必须一致。

`cooler zoomify` 通常从最细的分辨率逐层合并得到较粗的分辨率，因此较粗的层是较细层的整数倍。这只是生成方式的惯例；OmicsWeave 转换时直接读取每一层已有的像素，不依赖层与层之间的关系。

与单个 cool 一样，各层可能各自带有 `bins/weight` 等 balance 列，它们是分析结果，不影响原始 `count`。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 格式识别 | `formats.detect` 按 `/resolutions` 组识别为 `Mcool` |
| 区域查询 | 不直接支持，转换为 `.hic` 后用 `ContactReader` |
| 转换 | `.mcool` → `.hic`（v9），默认转换全部分辨率 |
| 选择分辨率 | `resolutions` / `--resolutions` 选择文件中已有的分辨率 |
| 每一层的要求 | 与 [cool](cool.md) 相同：layout v2 / v3、固定宽度 bin |
| genome 名称 | 来自 cooler 的 `genome-assembly` 属性，缺失时写 `unknown` |

不直接读取 mcool 的区域，不写 mcool，也不从 `.hic` 反向转换。

## 转换

```sh
# 转换全部已有分辨率
ow convert contacts.mcool contacts.hic -t 4 -m 2G

# 仅选择文件中已有的两个分辨率
ow convert contacts.mcool selected.hic --resolutions 5000,25000
```

```python
from omicsweave import convert
from omicsweave.query import ContactReader

convert.run("contacts.mcool", "selected.hic", resolutions=[5000, 25000])

with ContactReader("selected.hic") as contacts:
    print(contacts.resolutions())  # (25000, 5000)
```

- 选择不存在的分辨率会报错，错误信息列出文件中已有的分辨率；转换不会计算新的分辨率。
- 输出 `.hic` 中的分辨率按由粗到细排列，与 Juicer 工具一致。
- 每一层只使用上三角的原始 count，不使用 balance weight；square 存储的层同样只取上三角。
- 输出不带 expected 和归一化向量。

全分辨率的 mcool 转换时间和临时空间随最细分辨率的像素数增长。训练只需要一两种分辨率时，用 `--resolutions` 选择可以明显减少转换量和输出大小。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 不知道文件有哪些分辨率 | 直接转换一次并查看 `ContactReader.resolutions()`，或用 `cooler ls file.mcool` |
| 报分辨率不在文件中 | `--resolutions` 只能选已有分辨率，不能重新分箱 |
| 报各分辨率染色体不同 | 文件各层的染色体表不一致，需要先修正源文件 |
| 报 bin 宽度不固定 | 某层使用可变宽度 bin，当前不支持 |

## 相关格式

- [cool](cool.md)：单分辨率 cooler，mcool 的每一层都是一个 cool。
- [.hic](hic.md)：转换目标，同样在一个文件中保存多种分辨率。
- [4DN pairs](pairs.md)：分箱之前的接触对，也可以直接转换为指定分辨率的 `.hic`。

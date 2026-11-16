# cool

[格式总览](README.md) · [转换](../convert.md) · [.hic](hic.md) · [mcool](mcool.md) · [4DN pairs](pairs.md)

`.cool` 是 cooler 项目定义的**单分辨率接触矩阵**格式，建立在 HDF5 之上。它用几张列式表保存染色体、bin 和非零像素，是 Python Hi-C 生态（cooler、cooltools、HiCExplorer 等）最常用的存储方式。

OmicsWeave 不直接按窗口查询 `.cool`，而是把它转换为 [`.hic`](hic.md) 后再用 `ContactReader` 读取。转换只搬运原始计数，保留源文件的分辨率，不做归一化。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 一个分辨率的接触矩阵（稀疏像素） |
| 常见后缀 | `.cool` |
| 存储 | HDF5 文件，列式数据集，通常 gzip（deflate）压缩分块 |
| 坐标 | bin 表中的 `start` / `end` 是 0-based 半开碱基坐标 |
| 随机访问 | 格式本身支持（`indexes` 组）；OmicsWeave 用于转换，不直接查询 |
| 版本 | cooler layout v2 / v3 |
| 典型来源 | `cooler cload` / `cooler load`、distiller、4DN 数据门户 |
| OmicsWeave | 识别、转换为 `.hic` |

## 格式由来

cooler（Abdennur & Mirny, *Bioinformatics* 2020）把接触矩阵看成一张稀疏表：每个非零格子一行，记录两个 bin 编号和计数。表按列存放在 HDF5 中，可以用 NumPy / pandas 直接读取任意列。与 Juicer 的 [`.hic`](hic.md) 相比，cool 更像一个“数据库文件”：结构透明、易于用通用 HDF5 工具检查、方便附加自定义列（例如 balance weight）。

单个 `.cool` 只保存一种分辨率；多种分辨率合在一个文件里就是 [`.mcool`](mcool.md)。

## 文件结构

一个 `.cool` 文件的根组包含四个子组和若干属性：

```text
/
├── chroms/
│   ├── name       染色体名，例如 chr1、chr2 …
│   └── length     染色体长度
├── bins/
│   ├── chrom      每个 bin 所属染色体（编号或枚举）
│   ├── start      bin 起点
│   ├── end        bin 右端
│   └── weight     可选：balance weight（ICE 等校正向量）
├── pixels/
│   ├── bin1_id    行 bin 编号
│   ├── bin2_id    列 bin 编号
│   └── count      原始接触数
└── indexes/
    ├── chrom_offset   每条染色体第一个 bin 在 bins 表中的位置
    └── bin1_offset    每个行 bin 第一个像素在 pixels 表中的位置
属性：format-version、bin-type、bin-size、storage-mode、genome-assembly …
```

- **bins** 是基因组的分箱表。固定宽度分辨率下，每条染色体从 0 开始按 `bin-size` 切分，最后一个 bin 可能较短。
- **pixels** 只存非零格子，按 `(bin1_id, bin2_id)` 排序。`bin1_id`、`bin2_id` 是全基因组统一编号，可以跨染色体。
- **indexes** 让读取方无需扫描就能找到某条染色体的 bin 或某一行的像素，类似 CSR 稀疏矩阵的行指针。

### storage-mode：上三角或完整矩阵

| `storage-mode` | 存储内容 |
|---|---|
| `symmetric-upper`（默认） | 只存 `bin1_id <= bin2_id` 的格子，矩阵视为对称 |
| `square` | 存储完整矩阵的全部非零格子，可以表示非对称矩阵 |

### balance weight 与原始计数

`pixels/count` 是原始接触数。`cooler balance` 计算的 ICE 校正向量写入 `bins/weight`，校正后的值是 `count × weight[bin1] × weight[bin2]`。weight 是分析结果，矩阵本身仍是原始计数。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 格式识别 | `formats.detect` 识别为 `Cool`（单分辨率）或 [`Mcool`](mcool.md) |
| 区域查询 | 不直接支持，转换为 `.hic` 后用 `ContactReader` |
| 转换 | `.cool` → `.hic`（v9），保留源分辨率 |
| layout | cooler format-version 2 / 3 |
| bin 类型 | 仅固定宽度 bin；可变宽度 bin（例如按酶切片段分箱）会报不支持 |
| storage-mode | `symmetric-upper` 和 `square`；square 只取上三角 |
| genome 名称 | 来自 `genome-assembly` 属性，缺失时写 `unknown` |

转换时不读取 `bins/weight` 或其他 balance 列，输出 `.hic` 中没有 expected 和归一化向量。当前不提供 `.hic` → cool 的反向转换，也不写新的 cooler 文件。

## 转换与查询

```sh
ow convert contacts.cool contacts.hic -t 4 -m 2G
```

```python
from omicsweave import convert
from omicsweave.query import ContactReader
from omicsweave.types import ContactWindow

report = convert.run("contacts.cool", "contacts.hic", threads=4, memory="2G")
print(report.records, report.seconds)

with ContactReader("contacts.hic") as contacts:
    print(contacts.resolutions())  # 与源 cool 的 bin-size 相同
    window = ContactWindow(("chr1", 0, 1_000_000), ("chr1", 0, 1_000_000), contacts.resolutions()[0])
    matrix = contacts.matrix(window)
```

转换时像素并行读取（直接读 HDF5 的压缩块并在多个线程中解压），放入目标 `.hic` 的块中，在内存或磁盘临时文件里排序后逐块写出。`threads` / `-t` 控制线程数，`memory` / `-m` 控制缓冲预算，超出时使用临时文件。

### 转换后的含义

- 输出分辨率就是源 cool 的分辨率，不能用 `resolutions` 重新分箱；该选项只用于从 [mcool](mcool.md) 中选择已有分辨率，对 `.cool` 传入会报错。
- 计数是 `pixels/count` 的原始值，不乘 balance weight。
- `square` 模式只取上三角（`bin1_id <= bin2_id`），不会把上下三角相加；`.hic` 对每个格子只存一次。
- 查询 `.hic` 时，对称部分由 `ContactReader` 补全；对 `symmetric-upper` 文件，读回的矩阵就是 cooler 以 `balance=False` 读出的稠密矩阵。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 报 bin 宽度不固定 | 文件使用可变宽度 bin；当前只支持固定宽度分辨率 |
| 打开 `.mcool` 时报错 | 多分辨率文件使用 [mcool](mcool.md) 的转换方式 |
| 转换后的值比 cooler 显示的小 | cooler 显示可能是 balanced 值；OmicsWeave 保留原始 count |
| 想转换时改变分辨率 | 先用 cooler 生成所需分辨率（例如 `cooler zoomify`），再转换 |
| 需要 balance 后的矩阵 | 读出原始矩阵后，用 cooler 计算的 weight 自行相乘 |

## 相关格式

- [mcool](mcool.md)：同一矩阵的多分辨率 cooler 文件。
- [.hic](hic.md)：OmicsWeave 查询接触矩阵使用的格式。
- [4DN pairs](pairs.md)：分箱之前的接触对文本，cooler 的常见输入。

# 安装与快速开始

[文档首页](README.md) · [格式与索引](formats.md) · [区域查询](queries.md)

## 从源码安装

需要 Python 3.12+、Rust 1.89+、C 编译器和 CMake 3.26+。HDF5 从源码静态构建，因此初次安装包含 native 编译。

```sh
uv pip install 'git+https://github.com/CladeTeam/OmicsWeave'
```

已有仓库时，在根目录执行：

```sh
uv sync --locked
uv run ow --help
```

安装提供 `omicsweave` Python 包和 `ow` 命令。本文以源码安装为准；发布 wheel 的构建与安装验证见[开发与发布](development.md)。PyTorch、pyarrow 和基准对照库按用途另装，不是普通查询的必需依赖。

## 可直接运行的小例子

下面自行写一个小型 bgzip FASTA，writer 同时生成索引，不需要下载数据：

```python
from omicsweave import Region
from omicsweave.formats import FastaWriter
from omicsweave.query import SequenceReader

with FastaWriter("example.fa.gz") as writer:
    writer.add("chr1", "AACCGGTT")

with SequenceReader("example.fa.gz") as genome:
    key = Region("chr1", 1, 5)
    result = genome.fetch(key)
    print(result.text)                             # ACCG
    print(result.bases.dtype, result.bases.shape)   # uint8 (4,)
    print(genome.fetch(key.with_strand("-")).text)  # CGGT
```

`Region("chr1", 1, 5)` 表示 `[1, 5)`，长度为 4；返回值是 ASCII 碱基的 NumPy 数组，不是 one-hot 编码。默认结果在后续查询及 reader 关闭后仍然有效。

## 读取已有参考和信号

```python
from omicsweave import Region
from omicsweave.formats import BgzipFasta
from omicsweave.query import SequenceReader, TrackReader

# 对没有索引的现有文件执行一次；会遍历文件。
BgzipFasta("genome.fa.gz").build_index()

with SequenceReader("genome.fa.gz") as genome, TrackReader("atac.bw") as atac:
    key = Region("chr1", 1_000_000, 1_131_072)
    bases = genome.fetch(key).bases       # uint8 [131072]
    signal = atac.signal(key).values      # float32 [131072]，缺失为 NaN
    bins = atac.bins(key, 128, "mean")    # float32 [1, 1024]
```

Reader 不会在打开时建索引。普通 gzip FASTA 需先转换为 bgzip FASTA 或 2bit；VCF、BED、SAM 和 cooler 也需先[转换](conversion.md)才能使用对应的区域 reader。

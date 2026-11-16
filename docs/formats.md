# 格式与索引

[文档首页](README.md) · [区域查询](queries.md) · [格式转换](conversion.md)

## 区域读取与预处理

| 文件内容 | 区域 reader | 必需索引 / 参考 | 可接受的预处理源 |
|---|---|---|---|
| 普通 FASTA | `SequenceReader` | `.fai` | FASTA、gzip FASTA、BGZF FASTA |
| bgzip FASTA | `SequenceReader` | `.fai` + `.gzi` | 同上 |
| 2bit | `SequenceReader` | 内置目录，无外部索引 | FASTA 的各种存储 |
| bigWig | `TrackReader` | 内置索引 | bedGraph、wig、bedMethyl，含 gzip |
| bigBed | `IntervalReader` | 内置索引 | BED、bedMethyl，含 gzip |
| hic v8 / v9 | `ContactReader` | 内置索引 | cool / mcool、4DN pairs、Juicer short text |
| BCF（BGZF） | `VariantReader`、`ConsensusReader` | `.csi`；consensus 还需参考 reader | VCF，含 gzip / BGZF |
| BAM | `AlignmentReader` | `.bai` 或 `.csi` | SAM，含 gzip / BGZF；CRAM；BAM 排序 |
| CRAM 3.0 / 3.1 | `AlignmentReader` | `.crai`；按字段和文件模式可能需要参考 | 读取或转 BAM，不提供 CRAM writer |

`formats` 中的类描述文件，不代表每类都有区域 reader 或 writer。BED、VCF 文本、普通 gzip FASTA、SAM、bedGraph、wig、bedMethyl、cool / mcool 和 pairs 不直接用于区域查询。GTF / GFF3 由 [Annotation](annotation.md) 完整加载；[定量分析](quantification.md) 顺序读取坐标排序 BAM，不要求 BAM 区域索引。

格式按内容识别，扩展名用于转换目标推断。例如：

```python
from omicsweave import formats, open

file = formats.detect("sample.data")
print(file.format, file.path)
with open(file) as reader:
    print(reader.chroms)
```

`open` 选择单文件 reader；consensus 需显式创建 `ConsensusReader` 并传入参考。

## 显式构建索引

```python
from omicsweave.formats import Bam, BgzipFasta, Cram, Fasta

Fasta("genome.fa").build_index()             # .fai
BgzipFasta("genome.fa.gz").build_index()     # .fai + .gzi
Bam("sorted.bam").build_index()             # 默认 auto：BAI，过长参考用 CSI
Cram("reads.cram").build_index()            # .crai
Bcf("calls.bcf").build_index()              # .csi，与 bcftools index 相同
```

BCF 也可在转换时用 `ow convert variants.vcf.gz variants.bcf --index` 同时生成 CSI；已有的有序 BGZF BCF 用 `Bcf.build_index(min_shift=14)` 建索引（结果与 `bcftools index` 相同，乱序记录报错）。建索引不会排序输入；BAM 排序用 `ow convert in.bam sorted.bam --sort --index`。

| Reader | 默认查找次序 |
|---|---|
| FASTA | `<source>.fai`，BGZF 另需 `<source>.gzi` |
| BCF | `<source>.csi` |
| BAM | `<source>.csi`、去掉 `.bam` 的 `.csi`、`<source>.bai`、去掉 `.bam` 的 `.bai` |
| CRAM | `<source>.crai`、去掉 `.cram` 的 `.crai` |

索引不在旁边时可传格式对象，例如 `SequenceReader("genome.fa.gz", Fai("index/genome.fai"), Gzi("index/genome.gzi"))`，BAM / CRAM 的 `index=` 同样可显式指定。缺失必需索引会报 `IoError`。

## 写入

```python
import numpy as np
from omicsweave.formats import BigBedSchema, BigBedWriter, BigWigWriter, FastaWriter

chroms = {"chr1": 1000}
starts = np.array([100, 300], dtype=np.uint32)
ends = np.array([200, 400], dtype=np.uint32)
with FastaWriter("genome.fa.gz") as writer:
    writer.add("chr1", "ACGT" * 250)
with BigWigWriter("signal.bw", chroms) as writer:
    writer.add("chr1", starts, ends, np.array([1, 2], dtype=np.float32))
with BigBedWriter("peaks.bb", chroms, BigBedSchema(field_count=6)) as writer:
    writer.add("chr1", starts, ends, ["peak1\t0\t+", "peak2\t0\t-"])
```

Writer 在正常退出 `with` 时完成写入；异常退出不调用正常 finish。`close()` 显式完成并返回伴随文件，重复关闭返回空 tuple。Writer 是进程内对象，应在使用它的进程中新建。

`HicWriter` 接受已经分 bin 的像素，`HicBuilder` 从接触对构建各分辨率。更完整的文件级写入可使用[转换入口](conversion.md)。

## Sidecar

bigWig、bigBed 和 hic 的 `.sidecar` 保存块的尺寸和布局，供外部加载器规划内存。通过文件对象 `build_sidecar()` 或 writer 的 `sidecar=True` 创建；普通区域查询无需 sidecar，`BlockCache` 也不要求它。

Sidecar 按源文件长度和抽样内容指纹校验。源文件更新后应重建，不应把抽样指纹视为全文件加密校验。二进制布局见[Rust 接口](rust-api.md#sidecar)。

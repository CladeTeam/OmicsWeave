# 格式转换

[按区域查询](query.md) · [定量分析](quant.md)

`convert` 把现有文件转换为 reader 使用的查询格式，包括基因组数据和 FFM motif 集合，也支持 BAM 坐标排序。命令形式为 `ow convert 输入 输出`。

```sh
ow convert genome.fa genome.fa.gz -t 4 -m 2G
ow convert variants.vcf.gz variants.bcf --index
ow convert reads.sam reads.bam
ow convert reads.bam sorted.bam --sort --index
```

转换源默认按内容识别，目标默认按扩展名选择。上面的 FASTA 输出是 bgzip，并自动生成 `.fai` 和 `.gzi`；`--index` 为 BCF / BAM 额外建索引。示例中的输入、输出路径需要替换为自己的文件。

## 支持的转换

下表中的文本输入均支持纯文本、普通 gzip 和 BGZF 压缩。转换方向仅限表中列出的组合。

| 输入 | 输出 | 索引或伴随文件 | 说明 |
|---|---|---|---|
| [FASTA](formats/fasta.md) | bgzip FASTA | 自动 `.fai`、`.gzi` | 包括重新压缩已有 gzip / BGZF FASTA |
| FASTA | [2bit](formats/2bit.md) | 内置目录 | 非 ACGT 碱基转为 N，保留软屏蔽信息 |
| [VCF](formats/vcf.md) | [BCF](formats/bcf.md) | `index=True` 时 `.csi` | 不自动排序 |
| [SAM](formats/sam.md) | [BAM](formats/bam.md) | `index=True` 时 `.bai` / `.csi` | 不自动排序 |
| BAM | 按坐标排序的 BAM | 可选 `.bai` / `.csi` | 必须指定 `sort=True` |
| [CRAM](formats/cram.md) | BAM | 可选 `.bai` / `.csi` | 通常需要压缩时的参考，不自动排序 |
| [BED](formats/bed.md)、[bedGraph](formats/bedgraph.md) | [bigBed](formats/bigbed.md) | 内置索引，可选 `.sidecar` | 保留坐标及附加字段 |
| bedGraph、[wig](formats/wig.md) | [bigWig](formats/bigwig.md) | 内置索引，可选 `.sidecar` | 存储数值信号 |
| [bedMethyl](formats/bedmethyl.md) | bigBed | 内置索引，可选 `.sidecar` | 保留全部列及方言 schema |
| bedMethyl | bigWig | 内置索引，可选 `.sidecar` | 选择一个数值列，可按修饰、链筛选 |
| [motif 文本](formats/motif-texts.md)（MEME、JASPAR、TRANSFAC、HOMER、CIS-BP、HOCOMOCO、Cluster-Buster、BaMM、ProBound、GLAM2）、[FFM](formats/ffm.md) | FFM | 自带目录 | 可多个来源文件；按各格式规则保留计数、概率或权重模型 |
| [MAF](formats/maf.md) | [bigMaf](formats/bigmaf.md)（bigBed） | 内置索引，可选 `.sidecar` | 每个比对块一条 bed3+1 记录，按 `mafToBigMaf` 的文本布局写入 |
| [cool](formats/cool.md)、[mcool](formats/mcool.md) | [.hic](formats/hic.md) | 内置索引，可选 `.sidecar` | 原始 count，不使用 balance weight |
| [4DN pairs](formats/pairs.md)、[Juicer short](formats/juicer-short.md)、[HiC-Pro valid pairs](formats/hicpro-pairs.md) | `.hic` | 内置索引，可选 `.sidecar` | 按分辨率汇总接触 |

BAM → bigWig、BAM → 计数矩阵以及 BAM → Hi-C 接触属于 [quant](quant.md)。当前 convert 不提供 BCF → VCF、BAM → SAM、BAM → CRAM、2bit → FASTA 或 hic → cool 等反向导出。

## Python 调用与返回值

```python
from omicsweave import convert

report = convert.run(
    "cohort.vcf.gz", "cohort.bcf",
    index=True, threads=4, memory="2G",
)
print(report.target)      # 主输出格式对象
print(report.companions)  # 索引、sidecar 等伴随文件
print(report.files)       # 主文件在前的全部输出
print(report.records, report.seconds)
print(report.warnings)
```

`report` 还提供 `source`、`bytes_read`、`spilled`。`records` 的单位随转换而变，可能是序列、变异、比对、区间、像素或 motif；`spilled` 是写入临时 spill 文件的字节数。跳过未知染色体等情况记录为 `warnings`，Python 同时发出 `ConvertWarning`，CLI 打印到 stderr。

## 格式识别和显式指定

```python
from omicsweave.formats import Bed, BigBed, detect

print(detect("input.dat"))
convert.run(Bed("input.dat"), BigBed("output.dat"), chrom_sizes="genome.chrom.sizes")
```

```sh
ow convert input.dat output.dat --from bed --to bigbed --chrom-sizes genome.chrom.sizes
```

创建格式对象只是声明路径和预期格式，不会创建文件。四列且最后一列为数字的 BED 可能被自动识别为 bedGraph，需要保留区间语义时显式指定 `Bed` / `--from bed`。

`convert.routes()` 返回支持的源、目标格式类对；`ow convert --help` 列出格式名与全部选项。不适用于当前转换的选项会报错。

## FASTA 与索引

```sh
ow convert genome.fa genome.fa.gz
ow convert genome.fa.gz genome.2bit
```

`.fa.gz` 目标使用 bgzip，写出后可直接交给 `SequenceReader`。普通 gzip FASTA 作为输入可以转换，但不能直接做区域查询。FASTA → 2bit 会将 ACGT 以外的字符存为 N，因此需要保留完整 IUPAC 字符时应使用 FASTA。

只想给已有 FASTA 建索引时，使用格式对象的 `build_index()`，不必转换；方法见 [FASTA](formats/fasta.md)。

## VCF、SAM、CRAM 和 BAM 排序

```sh
ow convert cohort.vcf.gz cohort.bcf --index
ow convert coordinate-sorted.sam reads.bam --index
ow convert unsorted.bam sorted.bam --sort --index
ow convert reads.cram reads.bam --reference genome.fa.gz --index
```

VCF → BCF、SAM → BAM 和 CRAM → BAM 保留输入记录顺序。建索引要求输入已按坐标排序；`--sort` 仅用于 BAM → BAM。SAM 无序时，先转 BAM，再单独排序。当前没有 VCF 排序路由。

| Python 选项 | CLI | 用法 |
|---|---|---|
| `index=True` / `"auto"` | `--index` | BCF 写 CSI；BAM 在参考长度允许时写 BAI，否则写 CSI |
| `index="bai"` / `"csi"` | `--index bai` / `csi` | 为 BAM 显式选择索引 |
| `min_shift=N` | `--min-shift N` | 仅用于输出 BAM 且显式指定 `index="csi"` / `--index csi` 时，设置最小 bin 大小 `2**N`，默认 14；`auto` 即使最终选 CSI 也不接受。VCF → BCF 的 CSI 固定为 14，传入会报错（已有 BCF 可用 `Bcf(...).build_index(min_shift=N)`） |
| `add_contigs=True` | `--add-contigs` | 为 VCF header 缺少的 contig 补声明，多扫描一次输入 |
| `reference="genome.fa.gz"` | `--reference` | CRAM 解码使用的参考，需匹配原始压缩参考 |
| `verify_checksums=False` | `--no-verify` | CRAM 跳过块 CRC32 和参考 MD5 检查 |

VCF header 中未声明的 INFO / FORMAT 等名字应先修正，`add_contigs` 只处理 contig。SAM 的 RNAME 不在 `@SQ` 中时会报错。BAM 排序按坐标进行，不添加 `@PG`；已有排序 BAM 若只缺索引，直接用 `Bam(...).build_index()`，见 [BAM](formats/bam.md)。

## 区间和信号

```sh
ow convert peaks.bed peaks.bb --chrom-sizes genome.fa.gz.fai
ow convert signal.bedGraph signal.bw --chrom-sizes genome.chrom.sizes
ow convert signal.wig signal.bw --chrom-sizes genome.2bit
```

`chrom_sizes` / `--chrom-sizes` 接受两列染色体长度表、`.fai` 或 `.2bit`。不提供时，以每条染色体观测到的最大区间右端作为长度；需要使用完整参考长度时应显式提供。

这些路由会排序，不要求文本事先排序。bigBed 可以保存重叠区间；bigWig 输入区间不能重叠。bigWig 转换拒绝空区间和非有限信号值；超过声明染色体长度的区间会报错，不会自动裁剪。

wig 支持 `fixedStep`、`variableStep` 和 bedGraph 行，源坐标按各自语法解释。`track` / `browser` 行会跳过；fixedStep 的 span 大于 step 而造成重叠时会报错。

输出布局可用 `bbi=BbiOptions(...)` 设置 blocks、索引和 zoom；CLI 提供 `--block-size`、`--items-per-slot`、`--zoom-levels`。日常转换可使用默认值。

## bedMethyl

```sh
ow convert calls.bed calls.bb
ow convert calls.bed methylation.bw --mod m --value percent --strand +
```

支持 modkit 18 列和 ENCODE WGBS 11 列，两种方言不能混在同一输入中。转 bigBed 保留全部字段；转 bigWig 时：

- `value="percent"`（默认，也可写 `percent_modified`）：原始百分数，仍是 0–100。
- `value="coverage"`（也可写 `valid_coverage`）或 `"score"`：保存对应列。
- modkit 还支持 `nmod`、`ncanonical`、`nother_mod`、`ndelete`、`nfail`、`ndiff`、`nnocall`。
- 多个修饰代码并存时必须指定 `mod`；两条链的条目重叠时，需通过 `strand` 选择。`strand="."` 选择标记为合链的记录，不会替你合并正负链。

筛选后仍重叠会报错。转换不从计数重新计算百分比，也不把百分数换成 0–1 的 beta。字段说明见 [bedMethyl](formats/bedmethyl.md)。

## 多序列比对（MAF → bigMaf）

```sh
ow convert chr22.maf.gz multiz100way.chr22.bigmaf
ow convert multiz.maf.gz multiz.bb --ref-species hg38 --chrom-sizes hg38.chrom.sizes
```

```python
convert.run("multiz.maf.gz", "multiz.bigmaf", ref_species="hg38", threads=8)
```

目标是 bigMaf，它就是 bigBed（扩展名 `.bigmaf`、`.bb`、`.bigbed` 均可）：每个 `a` 比对块一条记录，坐标取参考物种那一行，第 4 列 `mafBlock` 保存整个块。转换采用 UCSC `mafToBigMaf` 的文本布局：`a score=` 保留六位小数，各列按块内最宽的值对齐，每行以 `;` 结尾；autoSql 使用 `bigMaf.as` 的定义。

- `ref_species` / `--ref-species`：参考物种。不指定时取第一个块首个 component 的 `src` 在第一个点之前的部分；参考名包含一个点时需显式指定，例如 `hg19.1`，转换器可匹配到第二个点。每个块都必须以参考物种的 component 开头，且在 `+` 链上，否则报错并给出行号。
- 染色体长度默认取参考行的 `srcSize`；同一染色体前后的 `srcSize` 不一致会报错。给出 `chrom_sizes` 时逐条核对，不一致同样报错。
- 输入不要求有序，按染色体（首次出现的顺序）和起点排序；块在参考上相互重叠会报错。
- 注释行、`##maf` 头和非 `a` 段落不进入输出；`##maf` 头必须存在并带 `version=`。

转换只改写格式，不做投影或过滤。按区域读取为 `[物种, 长度]` 矩阵见 [按区域查询](query.md) 的 `MsaReader`，格式细节见 [MAF](formats/maf.md) 与 [bigMaf](formats/bigmaf.md)。

## Motif 集合（→ FFM）

```sh
ow convert JASPAR2026_CORE_vertebrates_redundant_pfms_jaspar.txt jaspar.ffm
ow convert cisbp/pwms/*.txt cisbp.ffm
ow convert hocomoco/pfm/*.pfm more.pcm all.ffm --from pcm
```

```python
from pathlib import Path

from omicsweave import convert

convert.run(sorted(Path("cisbp/pwms").glob("*.txt")), "cisbp.ffm")
convert.run(["a.meme", "b.meme"], "all.ffm")
```

来源可以是一个文件或非空的文件序列，合并写成一个 FFM：按给出的顺序，相同路径重复传入时只处理一次；所有文件必须是同一格式、同一字母表，motif ID 必须全局唯一。CLI 中 shell 展开的 glob 就是多个来源，最后一个参数是目标文件；Python 中传路径列表（例如 `sorted(Path(...).glob(...))`）。OmicsWeave 不展开通配符。多源转换的 `plan.source`、`report.source` 表示第一个来源，完整输入列表由调用方保留。

计数、概率按原值保存；对数权重及能量系数转为 `exp(值)` 的 potential。每个 motif 保留 ID 和一个注释键 `name`，FFM → FFM 同样只带出 `name` 注释。没有行的 CIS-BP 空壳文件跳过并警告；大量跳过会合并为汇总警告。格式识别在 4 × 4 方阵等无法区分的情况下报错，需要 `--from` 指定；格式识别不会借助输入扩展名。各格式的读取规则见 [Motif 文本格式](formats/motif-texts.md)。

## 接触数据

```sh
ow convert sample.cool sample.hic
ow convert sample.mcool sample.hic --resolutions 5000,25000
ow convert sample.pairs.gz sample.hic --resolutions 5000,25000
ow convert sample.allValidPairs sample.hic --chrom-sizes genome.sizes
ow convert sample.short sample.hic --from juicer-short \
  --chrom-sizes genome.chrom.sizes --resolutions 5000,25000
ow convert sample.allValidPairs sample.hic --chrom-sizes genome.chrom.sizes \
  --resolutions 5000,25000
```

cool 读取已有分辨率，mcool 默认转换全部分辨率，也可选择其中一部分；不额外计算新的分辨率。只使用上三角原始 count，square 存储同样如此，不应用 balance weight。

pairs 的染色体信息默认来自 `#chromsize`；`chrom_sizes` 可覆盖。Juicer short 和 HiC-Pro validPairs / allValidPairs 没有该 header，必须提供染色体长度。HiC-Pro 接受 12 列和带可选 allele 的 13 列，允许行尾 tab；方向列为 `+` / `-`，每行计一次。HiC-Pro 的位置是 1-based 读段 5′ 端；转换按原值分箱，不减一，与 `hicpro2juicebox.sh` 交给 Juicer `pre` 的位置一致，详见 [HiC-Pro valid pairs](formats/hicpro-pairs.md)。

pairs → hic 的 `resolutions` 默认是 2.5M、1M、500k、250k、100k、50k、25k、10k、5k。指定的较粗分辨率必须是最细分辨率的整数倍。只写 BP 分辨率，结果为原始计数，不包含 expected、归一化向量或 Juicer 的 `All` 染色体。

每行按一次接触计数，short-with-score 使用该行 score；不会自动做 MAPQ、pair type、重复或同一酶切片段过滤。位置按输入值直接分箱，不减一；超出染色体长度的位置先截断为染色体长度再分箱（与 Juicer `pre` 相同），因此染色体长度恰为分辨率整数倍时，这些接触落在起点等于染色体长度的那个 bin。未知染色体行跳过并发出 warning。需要改变这些规则时，先清理输入，再转换。

## 资源预算、试运行和输出

```sh
ow convert unsorted.bam sorted.bam --sort --index \
  -t 4 -m 2G --tmp-dir scratch --dry-run
```

```python
plan = convert.plan(
    "unsorted.bam", "sorted.bam",
    sort=True, index=True, threads=4, memory="2G", tmp_dir="scratch",
)
print(plan.source, plan.target, plan.budget, plan.memory)
report = plan.run()
```

`plan` / `--dry-run` 检查转换组合、选项和内存规划，不生成转换输出，不代替完整输入校验。自动识别输入仍可能读取文件头；cool / mcool 规划还会打开选中的每个分辨率组，读取并校验完整 bins 表，以估算固定内存，因此不保证只读少量头信息。默认线程数是进程可用的 CPU 核心数，默认缓冲预算为 4 GiB；`memory` 接受字节数或 `"2G"` 这样的大小字符串。

支持 spill 的排序路由在缓冲不足时写临时文件，默认放在输出旁边，也可用 `tmp_dir` / `--tmp-dir` 指定。预算不包括所有分配器开销和线程栈，不是 RSS 硬上限。输出内容不随线程数或是否 spill 改变。

Motif → FFM 当前按文件整体解析，每次保留一个来源的模型；FFM 来源先读成 resident 集合，各来源的压缩输出块保留到最终写出。这条路由不做磁盘 spill，`threads` 也不将文件解析变成并行流水线。计划中的固定内存是源文件总字节数的四倍估算，压缩输入可能明显低估展开后的模型；运行时会拒绝超过 `memory` 的展开文本，并给 FFM 解码设置分配上限，但这些检查不等于整次转换的 RSS 硬限制。大集合应按实际内存占用选择预算。

`level` / `--level` 设置压缩等级：bgzip FASTA、BCF、BAM、bigBed、bigWig 为 0–12；hic 为 0–9，大于 9 时按 9 处理。FASTA → 2bit 和 motif → FFM 不接受 `level`。`sidecar=True` / `--sidecar` 可为 bigBed、bigWig、hic 额外写 [.sidecar](formats/sidecar.md)，普通区域查询不需要它。

所有结果先写临时文件，成功后发布主文件和伴随文件。多文件发布不是文件系统事务，应等 `run()` 成功返回后再使用整组文件。CLI 成功时将输出路径打印到 stdout，进度和警告打印到 stderr。

## 示例：准备一套训练数据

下面把一组常见的公开数据整理成可以直接被 reader 查询的文件。每一步都只需执行一次，之后训练时只读转换结果。

```sh
mkdir -p ref tracks peaks align variants contacts reads   # 输出目录需要事先存在

# 参考基因组：普通 gzip → bgzip，自动生成 .fai 和 .gzi
ow convert hg38.fa.gz ref/hg38.bgz.fa.gz -t 8

# 信号轨道：bedGraph / wig → bigWig，用参考的 .fai 声明染色体长度
ow convert atac.bedGraph.gz tracks/atac.bw --chrom-sizes ref/hg38.bgz.fa.gz.fai
ow convert phyloP.wig.gz tracks/phyloP.bw --chrom-sizes ref/hg38.bgz.fa.gz.fai

# peaks：BED → bigBed，保留名称、分值等附加列
ow convert peaks.narrowPeak.gz peaks/atac.bb --chrom-sizes ref/hg38.bgz.fa.gz.fai

# 甲基化：modkit pileup → 正链 5mC 百分比 bigWig
ow convert pileup.bed.gz tracks/5mC.plus.bw --mod m --strand + --value percent

# 变异：VCF → BCF + CSI，可用于基因型查询和 consensus
ow convert cohort.vcf.gz variants/cohort.bcf --index -t 8

# 多序列比对：MAF → bigMaf，用于保守性 / MSA 模型
ow convert multiz100way.maf.gz align/multiz100way.bigmaf --ref-species hg38 -t 8

# 接触：mcool → .hic，仅保留需要的两个分辨率
ow convert sample.mcool contacts/sample.hic --resolutions 5000,25000 -t 8 -m 8G

# 比对：未排序 BAM → 排序 + 索引
ow convert aligned.bam reads/sorted.bam --sort --index -t 8 -m 8G --tmp-dir scratch
```

之后：

```python
from omicsweave.query import ContactReader, IntervalReader, SequenceReader, TrackReader

genome = SequenceReader("ref/hg38.bgz.fa.gz")
atac = TrackReader("tracks/atac.bw")
peaks = IntervalReader("peaks/atac.bb")
contacts = ContactReader("contacts/sample.hic")
```

大批量转换适合写成脚本，用 `convert.plan(...)` 先检查每一项的组合和选项，再逐项 `plan.run()`。转换阶段失败会清理临时输出；发布阶段逐个重命名文件，若中途失败，可能已发布部分文件，应检查主文件和伴随文件是否齐全。

## CLI 选项参考

| 选项 | 适用转换 | 作用 |
|---|---|---|
| `--from FORMAT` / `--to FORMAT` | 全部 | 显式指定源 / 目标格式，覆盖内容识别和扩展名推断 |
| `-t N`, `--threads N` | 全部 | 并行阶段的线程数，默认是进程可用的全部核心 |
| `-m SIZE`, `--memory SIZE` | 全部 | 缓冲预算，如 `512M`、`8G`，默认 `4G` |
| `--tmp-dir DIR` | 全部 | spill 文件位置，默认在目标文件旁 |
| `-l L`, `--level L` | 有压缩等级选项的路由 | bgzip、BCF、BAM、bigBed、bigWig 为 0–12；hic 为 0–9，较大值按 9 处理。2bit 和 FFM 目标不接受 |
| `-n`, `--dry-run` | 全部 | 检查组合、选项和预算并打印计划，不写文件 |
| `-q`, `--quiet` | 全部 | 只打印写出的文件 |
| `--index [auto\|bai\|csi]` | → BCF、→ BAM | 同时建索引；BCF 总是 CSI，BAM 的 `auto` 在参考长度允许时写 BAI |
| `--min-shift N` | → BAM，且 `--index csi` | CSI 最小 bin 为 `2^N`，默认 14；`--index auto` 不接受 |
| `--add-contigs` | VCF → BCF | 为 header 缺少的 contig 补声明，多扫描一遍输入 |
| `--sort` | BAM → BAM | 按坐标排序，与 `samtools sort` 相同 |
| `--reference FASTA` | CRAM → BAM | 压缩时使用的参考：FASTA + `.fai`、bgzip FASTA + `.fai` / `.gzi` 或 2bit |
| `--no-verify` | CRAM → BAM | 跳过块 CRC32 与 slice 参考 MD5 检查 |
| `--chrom-sizes FILE` | → bigBed、→ bigWig、pairs / short / HiC-Pro → hic | 染色体长度（长度表、`.fai` 或 `.2bit`）；MAF 默认取参考行的 `srcSize` 并据此核对；pairs 默认读 `#chromsize`，Juicer short 和 HiC-Pro 必须提供 |
| `--resolutions R,R,...` | mcool / pairs / short / HiC-Pro → hic | mcool：选择已有分辨率（默认全部）；接触对文本：分箱大小，默认 2.5M 到 5k |
| `--sidecar` | → bigBed、→ bigWig、→ hic | 同时写数据加载用的 `.sidecar` |
| `--block-size N` | → bigBed、→ bigWig | 索引节点的子节点数，默认 256 |
| `--items-per-slot N` | → bigBed、→ bigWig | 每个数据块的条目数，默认 bigBed 512、bigWig 1024 |
| `--zoom-levels N` | → bigBed、→ bigWig | 最多生成的 zoom 层数，默认 10 |
| `--value COLUMN` | bedMethyl → bigWig | 写入的数值列，默认 `percent` |
| `--mod CODE` | bedMethyl → bigWig | 修饰代码（第 4 列，如 `m`、`h`、`a`）；存在多种修饰时必须指定 |
| `--strand {+,-,.}` | bedMethyl → bigWig | 只取该链的条目；默认全部，链间重叠会报错 |
| `--ref-species NAME` | MAF → bigMaf | 参考物种（如 `hg38`），默认取第一个块第一行的物种 |

`ow convert --help` 末尾列出全部转换组合；Python 中用 `convert.routes()` 得到相同的（源格式类，目标格式类）列表。Python 关键字与 CLI 选项一一对应（`--chrom-sizes` → `chrom_sizes=`，`--no-verify` → `verify_checksums=False`，bigBed / bigWig 布局选项通过 `bbi=BbiOptions(...)` 传入）。

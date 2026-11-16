# BCF

[格式总览](README.md) · [VCF](vcf.md) · [CSI](csi.md) · [查询](../query.md) · [转换](../convert.md) · [定量](../quant.md)

BCF（Binary Call Format）是 [VCF](vcf.md) 的二进制形式：同样的 header、同样的变异记录和样本基因型，但字段预先编码成类型化的二进制，并以 BGZF 分块压缩。配合 [CSI](csi.md) 索引，可以只解码目标区域、只解码需要的字段。

在 OmicsWeave 中，BCF 是变异数据的标准查询格式：`VariantReader` 按区域返回基因型矩阵和 INFO / FORMAT 字段，`ConsensusReader` 把变异应用到参考序列上生成个体化序列，`quant` 用它屏蔽已知变异位点和统计等位计数。VCF 先[转换](../convert.md)为 BCF 再使用。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 与 VCF 相同：位点、等位、INFO、样本 FORMAT / GT |
| 常见后缀 | `.bcf` |
| 存储 | BGZF 压缩的二进制记录；header 为 VCF 文本 |
| 坐标 | 文件内 POS 为 0-based；OmicsWeave 结果为 0-based 半开 |
| 随机访问 | BGZF 压缩并按坐标排序时，通过 CSI 支持 |
| 外部索引 | `.csi`（BCF 不使用 BAI 或 tabix） |
| 规范 | hts-specs BCF2（2.1 / 2.2） |
| 典型来源 | `bcftools view -Ob`、`bcftools call -Ob`、`ow convert x.vcf.gz x.bcf` |
| OmicsWeave | 区域查询（基因型、字段、过滤表达式）、consensus、建 CSI、quant 辅助输入 |

## 格式由来

BCF 与 VCF 同属 hts-specs，由 htslib / bcftools 实现并广泛使用。早期的 BCF1 与 samtools mpileup 绑定，已不再使用；现行的是 BCF2，文件以 magic `BCF\2\1`（2.1）或 `BCF\2\2`（2.2）开头。

它解决的是文本 VCF 在大队列上的效率问题：一个位点有十万个样本时，文本需要切分十万个字段、逐个解析数字；BCF 中这些值是连续存放的定宽整数，可以整块读取，还可以跳过不需要的字段。

## 文件结构

```text
┌────────────────────────────────────┐
│ magic "BCF\2\2"                    │
│ header 长度 + VCF header 文本      │  与 VCF 的 ## 行和 #CHROM 行相同
├────────────────────────────────────┤
│ 记录 1：共享部分 │ 样本部分        │
│ 记录 2：共享部分 │ 样本部分        │  整个流按 BGZF 分块压缩
│ ...                                │
└────────────────────────────────────┘
```

### 字典：用编号代替名字

BCF 记录中不写字段名和染色体名，而写编号：

- **contig 字典**：按 header 中 `##contig` 的顺序编号，记录中的 CHROM 是这个编号。
- **字符串字典**：FILTER、INFO、FORMAT 共用一个字典。`PASS` 固定为 0，其余按声明顺序编号；header 行中的 `IDX=` 可显式指定编号。

所以 VCF 中用到的每个名字都必须在 header 中声明，否则无法编码；这也是 VCF → BCF 转换要求 header 完整的原因。

### 一条记录

**共享部分**保存位点信息：

| 字段 | 内容 |
|---|---|
| CHROM | contig 编号 |
| POS | 0-based 位置 |
| rlen | 参考跨度（REF 长度，或由 INFO/END 等给出） |
| QUAL | float32 |
| ID、REF / ALT | 类型化字符串 |
| FILTER | 字典编号向量 |
| INFO | 依次为“字典编号 + 类型化值” |

**样本部分**按 FORMAT 字段组织：每个字段先写字典编号和类型，再把**所有样本**的值连续写出。读取一个字段时可以直接跳过其他字段的整段数据，这正是只取 GT 很快的原因。

### 类型化值

每组值前有一个类型字节：低 4 位是类型，高 4 位是个数（15 表示个数另行给出）。

| 类型 | 编码 |
|---|---|
| 1 / 2 / 3 | int8 / int16 / int32，写入时取能放下的最小宽度 |
| 5 | float32 |
| 7 | char（字符串） |

### 缺失值与向量结束

定长存储需要两种特殊值：**缺失**（VCF 中的 `.`）和**向量结束**（这条记录或样本的值比最大宽度少，后面是填充）。

| 类型 | 缺失 | 向量结束 |
|---|---|---|
| int8 | `0x80` | `0x81` |
| int16 | `0x8000` | `0x8001` |
| int32 | `0x80000000` | `0x80000001` |
| float32 | 位模式 `0x7F800001`（NaN） | 位模式 `0x7F800002`（NaN） |

GT 单独编码：每个等位写为 `(编号 + 1) << 1 | 定相位`，`0` 表示缺失等位，向量结束值表示该样本的倍性低于本条记录的最大倍性。

OmicsWeave 把这些值统一映射为易用的常量，见下文“结果中的缺失值”。

## 关键概念

### 记录与区域的重叠

区域迭代使用 BCF 固定字段中保存的跨度 `[POS, POS + rlen)` 筛选记录，不为判断重叠重新读取 INFO/END 等字段。OmicsWeave 从 VCF 转换时，会根据 REF 长度及 INFO/END、INFO/SVLEN、FORMAT/LEN 等计算并保存 rlen。因此对正常编码了跨度的记录：

- 删除 `ACGT → A` 覆盖 4 个参考碱基，与其中任一碱基重叠的区域都会返回它。
- `<DEL>` 加 `END=200` 的结构变异，与其跨度内任意位置重叠的区域都会返回它。

### 倍性不一致

同一区域中可能同时有二倍体和单倍体基因型（例如 PAR 区附近、`0` 与 `0/1` 混合）。OmicsWeave 用区域内的最大倍性确定矩阵的第三维，较短的样本补 `GT_ABSENT`。

### 过滤表达式

bcftools 的 `-i` / `-e` 表达式（如 `QUAL>30 && INFO/DP>10`、`GT="het"`）是筛选变异的通用语言。OmicsWeave 用 Rust 实现 bcftools 1.23.1 的 `filter.c` 表达式语义（不含 Perl 过滤），不需要安装 bcftools。对参考实现具有已定义语义的输入进行结果对照；触发未初始化缓冲读取的输入会安全返回 `UnsupportedError`。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 格式识别 | 按解压后的 magic `BCF\2` 识别为 `Bcf`，支持 2.1 与 2.2 |
| 基因型矩阵 | `VariantReader.genotypes`：只解码固定字段和 GT |
| 字段读取 | `VariantReader.fetch`：QUAL、FILTER、INFO 与非 GT、非 Flag 的 FORMAT，一次扫描；GT 始终在 `records.genotypes` 中 |
| 样本选择 | `samples=` 按名字选择并按请求顺序返回；默认全部样本 |
| 过滤表达式 | `include` / `exclude`，bcftools 语法，二者互斥 |
| 个体化序列 | `ConsensusReader`：`bcftools consensus` 1.23 的 Rust 实现，在内存中生成 |
| 建索引 | `Bcf.build_index()` 写 `<bcf>.csi`，等同 `bcftools index` |
| 索引查找 | 默认 `<bcf>.csi`，再试去扩展名后的 `.csi`，之后同样尝试 `.bai`；找到 BAI 后仍因类型不符报错。可用 `index=` 指定 CSI |
| 缓存与并发 | 可接 `BlockCache`；解码释放 GIL；可 `clone()`、pickle |
| 生成 | VCF（纯文本 / gzip / BGZF）→ BCF；外部工具写出的 BGZF BCF 同样可用 |
| quant | `Task.variants` 屏蔽已知位点；`Alleles` 统计二等位 SNV 的 ref / alt / other |

限制：

- 区域查询只支持 BGZF 压缩、按坐标排序且有 CSI 的 BCF；未压缩 BCF（`bcftools -Ou`）不能建索引。
- reader 打开时不会自动建索引；没有任何候选索引时报 `IoError` 并提示 `Bcf.build_index()`；若只找到合法 BAI，则报 `InvalidError`：`a BCF is indexed by a CSI`。
- `VariantReader` 和 `ConsensusReader` 只接受连续 `Region`，不接受 `Spliced`。
- 当前不提供 BCF → VCF 导出，也没有从内存写 BCF 的 writer。

## 准备

```sh
ow convert cohort.vcf.gz cohort.bcf --index
```

已有排序的 BGZF BCF、只缺索引时：

```python
from omicsweave.formats import Bcf

csi = Bcf("cohort.bcf").build_index()  # Csi('cohort.bcf.csi')
csi = Bcf("cohort.bcf").build_index(min_shift=12, threads=4)
```

`min_shift` 默认 14（最小 bin 16 kb），索引层数按最长 contig 自动确定；`threads` 默认使用进程可用的全部核心解压。记录未排序时报 `InvalidError`。

## 查询基因型

```python
from omicsweave.query import VariantReader

with VariantReader("cohort.bcf") as variants:
    print(variants.chroms)   # contig → 长度；header 未给长度时为 None
    print(variants.samples)  # header 中的样本顺序

    gt = variants.genotypes(("chr1", 0, 100_000), samples=["s1", "s2"])
    print(gt.start, gt.end)                  # int64，0-based 半开，按 rlen
    print(gt.ids.to_list())                  # ['rs1', '.', ...]
    print(gt.ref.to_list(), gt.alt.to_list())  # ALT 多个时以逗号连接，无 ALT 为空字符串
    print(gt.gt.shape, gt.ploidy)            # (记录数, 2, 最大倍性)
    print(gt.phased)                         # bool [记录数, 样本数]
```

以上面 [VCF](vcf.md) 示例中的三条记录为例：

| 记录 | s1 | s2 | `gt[i]` | `phased[i]` |
|---|---|---|---|---|
| `A → G` | `0/1` | `1\|1` | `[[0, 1], [1, 1]]` | `[False, True]` |
| `ACGT → A` | `0\|1` | `./.` | `[[0, 1], [-1, -1]]` | `[True, False]` |
| `C → T,G` | `1/2` | `0` | `[[1, 2], [0, -2]]` | `[False, True]` |

- `gt` 是 `int8 [记录数, 样本数, 倍性]`：0 是 REF，1 是第一个 ALT。
- `GT_MISSING = -1` 是缺失等位（`.`），`GT_ABSENT = -2` 表示该样本没有这一位的等位（倍性较低），或记录没有 GT。
- `phased[i, j]` 在第一个之后的等位都已定相时为真；单倍体样本取其唯一等位的相位位（与 htslib 解码一致）。
- 没有记录的区域返回 `(0, 样本数, 0)` 形状的矩阵。

常见的下游转换：

```python
import numpy as np

g = gt.gt
called = (g >= 0)
alt_dosage = np.where(called, g > 0, 0).sum(axis=2)  # 每个样本携带的 ALT 等位数
missing = (g == -1).any(axis=2)                       # 有缺失等位的样本
```

## 读取字段与过滤

```python
with VariantReader("cohort.bcf") as variants:
    records = variants.fetch(
        ("chr1", 0, 100_000),
        samples=["s1"],
        info=["DP", "AF", "DB"],
        format=["AD"],
        include="QUAL>30 && INFO/DP>10",
    )
    print(records.qual)               # float32，'.' 为 NaN
    print(records.filters.to_list())  # 'PASS'、'LowQual'；多个以 ';' 连接，'.' 为空字符串
    print(records.info["AF"])         # float32 [记录数, width]
    print(records.info["DB"])         # bool [记录数]
    print(records.format["AD"])       # int32 [记录数, 样本数, width]
    print(records.genotypes.gt)       # 同 genotypes() 的结果
```

字段名必须在 header 中声明，否则报 `INFO/XX is not declared in the header`；同一个字段不能重复请求。

`format=["GT"]` 会报错：GT 已包含在 `records.genotypes` 中，应从这里读取。声明为 Flag 的 FORMAT 字段同样被拒绝；Flag 只支持 INFO 字段。

| 字段类型 | 返回 |
|---|---|
| INFO Flag | `bool [记录数]` |
| INFO Integer / Float | `int32` / `float32 [记录数, width]` |
| FORMAT Integer / Float | `int32` / `float32 [记录数, 样本数, width]` |
| 字符串 | `StringVector`：INFO 每条记录一项；FORMAT 按 `记录下标 * 样本数 + 样本下标` 排列 |

`width` 是本次查询结果中该字段的最大值个数。例如 `AF` 是 `Number=A`，区域中有一条三等位记录时 `width=2`，其他记录的第二列为填充。

### 结果中的缺失值

| 常量（`omicsweave.types`） | 值 | 含义 |
|---|---|---|
| `INT_MISSING` | `-2**31` | 整数缺失（`.`） |
| `INT_END` | `-2**31 + 1` | 填充：该记录 / 样本的值较少，或记录没有这个字段 |
| `FLOAT_MISSING` | 位模式 `0x7F800001` | 浮点缺失，读作 NaN |
| `FLOAT_END` | 位模式 `0x7F800002` | 浮点填充，读作 NaN |

浮点的两种情况都显示为 NaN，需要区分时查看位模式：

```python
from omicsweave.types import FLOAT_END, FLOAT_MISSING, INT_END, INT_MISSING

af = records.info["AF"]
bits = af.view("uint32")
is_missing = bits == FLOAT_MISSING
is_padding = bits == FLOAT_END

dp = records.info["DP"]
valid = (dp != INT_MISSING) & (dp != INT_END)
```

### 过滤表达式

`include` 保留表达式为真的记录（`bcftools -i`），`exclude` 去掉它们（`-e`），二者不能同时给出。支持 bcftools 的完整表达式语言，例如：

```python
variants.fetch(region, include='FILTER="PASS" && QUAL>=30')
variants.fetch(region, include="INFO/AF>0.01")
variants.fetch(region, exclude='TYPE="indel"')
```

`genotypes()` 不接受过滤表达式；需要过滤时使用 `fetch`，再从 `records.genotypes` 取基因型。

`fisher` 的两个 FORMAT 参数必须在 header 中声明为 `Number=R`，例如 `fisher(FMT/AD[:0-1], FMT/ADR[:0-1])`；仅有相同数量的值不能代替该声明。参考实现对不符合这个参数域的输入可能产生未定义结果，OmicsWeave 会拒绝读取未写入的结果值。

## 个体化序列：consensus

`ConsensusReader` 把 BCF 中的变异应用到参考序列的一个区域上，等同于对该区域运行 `bcftools consensus`，但在内存中直接返回数组：

```python
from omicsweave.query import ConsensusOptions, ConsensusReader, SequenceReader
from omicsweave.types import AllAlt, Haplotype, Samples

with SequenceReader("genome.fa.gz") as genome, ConsensusReader(
    "cohort.bcf", genome, options=ConsensusOptions(include="QUAL>20", chain=True)
) as cohort:
    region = ("chr1", 1000, 2000)
    hap1 = cohort.fetch(region, Haplotype("s1", "1"))
    print(len(hap1), hap1.text[:20], hap1.chain)

    batch = cohort.fetch_batch(region, [Haplotype("s1", "1"), Haplotype("s1", "2")])
    print(batch.offsets)          # 变长序列打包：第 i 条是 bases[offsets[i]:offsets[i+1]]
    print([s.text for s in batch])
```

| 等位选择 | 对应 bcftools | 含义 |
|---|---|---|
| `AllAlt()`（默认） | `-s -` | 每条记录取第一个 ALT，不看样本 |
| `Samples(["s1", ...])` | `-s s1,...` | 用这些样本的基因型生成 IUPAC 码 |
| `Haplotype("s1", "1")` | `-s s1 -H 1` | 取样本的某条单倍型 |
| `Record(i, alt=1)` | 无（相当于只含这条记录的文件加 `-s -`） | 只套用区域内第 `i` 条记录的这个 ALT，其他记录不套用 |
| `Variant(start, ref, alt)` | 同上 | 按 0-based 位置和 REF/ALT 找到那条记录，只套用它 |

后两种用于变异效应预测：每个变异单独生成 ALT 序列，窗口里的其他变异不会混进来；不应用 `iupac_codes`，但保留 mask、`absent`、标记和 chain 的行为。

`ConsensusReader.records(region)` 返回 `ConsensusRecords`，包含 `start` / `end`（`int64`，原始参考跨度）、`ref`、逗号连接的 `alt` 和布尔 `masked`。顺序是经过区域重叠规则、`include` / `exclude` 后的文件顺序，`Record.index` 从 0 开始，`Record.alt` 从 1 开始。mask 命中的记录仍占下标，只是不应用其 ALT；负链 key 不改变这个列表的顺序。`Variant` 忽略 REF / ALT 的 ASCII 大小写，重复匹配选文件中第一条；无匹配、记录下标或 ALT 编号越界均报错。实际示例见[查询文档](../query.md#个体化序列consensus)。

`Haplotype` 还支持 `R`、`A`、`LR`、`LA`、`SR`、`SA`、`I`、`<N>pIu` 等 bcftools `-H` 写法，写错时立即报 `ValueError`。

`ConsensusOptions` 对应 bcftools consensus 的选项：

| 选项 | bcftools | 作用 |
|---|---|---|
| `iupac_codes` | `-I` | 输出 IUPAC 码而不是单个等位 |
| `missing` | `-M` | 缺失基因型使用的字符 |
| `absent` | `-a` | 没有记录覆盖的位置使用的字符 |
| `mark_del` / `mark_ins` / `mark_snv` | `--mark-*` | 标记删除、插入、SNV |
| `mask` / `mask_with` | `-m` / `--mask-with` | 用 BED 区间屏蔽，默认以 `N` 替换 |
| `regions_overlap` | `--regions-overlap` | 区域边界上的记录如何取舍，默认 `"record"` |
| `chain` | `-c` | 同时返回参考 → consensus 的 chain 文本 |
| `include` / `exclude` | `-i` / `-e` | 只应用通过过滤的记录 |

要点：

- 区域使用参考坐标；插入和缺失会改变输出长度，批量结果按 `offsets` 打包，不能假定等长。
- 负链区域返回 consensus 的反向互补（IUPAC 码同样取互补，大小写保留）；chain 始终是正链。
- 每次 Python `fetch`、`records` 和非空 `fetch_batch` 都重新载入区域。`fetch_batch` 内部共享一次参考和记录解码，适合多样本、多单倍型和逐记录选择；跨调用的数据块复用由显式 `BlockCache` 提供。
- `reference` 可以是已打开的 `SequenceReader`（共享，不由 consensus 关闭），也可以是 FASTA、bgzip FASTA 或 2bit 路径（由 consensus 打开和关闭）。
- 返回数组由结果持有，之后查询、清空缓存或关闭 reader 不使它们失效。同一个 reader 在线程间串行使用其解码状态；`clone()` 共享文件、索引和缓存，使用独立状态并行解码。已打开的区域 reader 和缓存不能直接跨 fork 使用，pickle 则在接收进程重新打开且不携带缓存。

## 用于定量

```python
from omicsweave import quant

task = quant.Task(
    ["cells.bam"], "results",
    variants="known.bcf",                          # 屏蔽已知 SNP，避免误计为转换或甲基化
    products=[quant.Alleles.single_cell("sites.bcf")],
)
```

- `Task.variants`：代谢标记的转换计数和甲基化调用会跳过这些位点。
- `Alleles(sites)`：只统计 REF / ALT 都是单个 ACGT 的二等位 SNV，多等位和 indel 跳过并计入统计。
- `AlleleGenes(sample=...)`：使用该样本已定相的杂合位点，把基因计数分为 `hap1` / `hap2` / `both`。

两者都需要 BCF 配套 CSI。详见 [quant](../quant.md)。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| `no index of x.bcf found` | 缺少 CSI；`Bcf("x.bcf").build_index()` |
| 建索引报 `unsorted positions` | 记录未排序；先用外部工具排序再转换 |
| 找不到样本 | `samples` 中的名字不在 header；查看 `reader.samples` |
| `unknown chromosome` | 区域的 contig 不在 header；查看 `reader.chroms` |
| `gt` 中出现 `-2` | 该样本倍性低于区域最大倍性，或记录没有 GT |
| AF 等字段出现 NaN | 缺失或填充，用位模式区分 `FLOAT_MISSING` / `FLOAT_END` |
| consensus 长度与区域不同 | indel 改变了长度；用 `chain=True` 获取坐标映射 |

## 相关格式

- [VCF](vcf.md)：BCF 的文本形式，转换来源。
- [CSI](csi.md)：BCF 的区域索引。
- [FASTA](fasta.md)、[2bit](2bit.md)：consensus 的参考序列。
- [BAM](bam.md)：quant 等位计数与位点屏蔽的主输入。

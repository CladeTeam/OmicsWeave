# CRAM

[格式总览](README.md) · [查询](../query.md) · [转换](../convert.md) · [BAM](bam.md) · [CRAI](crai.md) · [SAM](sam.md)

CRAM 是**基于参考序列压缩**的比对格式。它保存的字段与 [SAM](sam.md) / [BAM](bam.md) 相同，但不直接存读段碱基，而是记录读段与参考基因组的差异；名称、质量、FLAG 等字段也按“数据系列”分别压缩。同样的数据，CRAM 通常明显小于 BAM，是大型测序项目（如 1000 Genomes、UK Biobank）的归档格式。

在 OmicsWeave 中，CRAM 可以通过 `AlignmentReader` 直接**按区域读取**，返回结构与 BAM 完全相同；也可以转换为 BAM，供 [quant](../quant.md) 使用。读取时只解压查询需要的数据系列：只要位置、FLAG、CIGAR 时，既不需要参考，也不解码碱基和质量。

## 速览

| 项目 | 说明 |
|---|---|
| 数据 | 读段比对记录，字段同 SAM |
| 常见后缀 | `.cram` |
| 存储 | container / slice / block 结构，按数据系列分别压缩 |
| 坐标 | OmicsWeave 结果为 0-based 半开 |
| 随机访问 | 坐标排序 + 索引后支持 |
| 外部索引 | [CRAI](crai.md)（`.cram.crai`） |
| 额外依赖 | 恢复碱基通常需要压缩时使用的参考序列 |
| 典型来源 | `samtools view -C`、GATK、大型测序项目的数据发布 |
| OmicsWeave | CRAM 3.0 / 3.1 区域查询、建 CRAI、CRAM → BAM |

## 格式由来

CRAM 最早由 EBI 提出（Fritz et al., *Genome Research* 2011），思路是：大多数读段与参考几乎相同，只记录差异就能大幅压缩。之后由 GA4GH 在 hts-specs 中标准化，htslib 和 samtools 原生支持。

- **CRAM 3.0** 是目前最广泛使用的版本，使用 gzip、bzip2、lzma 和 rANS 4x8 等编码。
- **CRAM 3.1** 保持同样的结构，新增 rANS Nx16、自适应算术编码、fqzcomp 质量值编码和读段名分词编码，压缩率更高。

## 文件结构

```text
┌──────────────────────────────────────────┐
│ 文件定义：magic "CRAM"、版本号、文件 ID   │
├──────────────────────────────────────────┤
│ header container：SAM header 文本         │
├──────────────────────────────────────────┤
│ data container                           │
│  ├─ container header（参考、起点、跨度、记录数、CRC32）│
│  ├─ compression header（各数据系列的编码方式）│
│  ├─ slice 1                              │
│  │   ├─ slice header（参考范围、参考 MD5） │
│  │   ├─ core block（位编码的数据系列）     │
│  │   └─ external blocks（各数据系列独立压缩）│
│  └─ slice 2 ...                          │
├──────────────────────────────────────────┤
│ data container ...                       │
├──────────────────────────────────────────┤
│ EOF container                            │
└──────────────────────────────────────────┘
```

| 层级 | 作用 |
|---|---|
| container | 一组 slice 的外壳，记录覆盖的参考区域，索引定位到这一层 |
| compression header | 说明本 container 中每个数据系列使用哪种编码、存在哪个 block |
| slice | 一批记录（htslib 默认每个 slice 约 10000 条），可独立解码 |
| block | 实际的压缩数据；每个 block 选择自己的压缩方法 |

### 数据系列

CRAM 不按记录逐条存储，而是把同一字段的所有值放在一起压缩，例如：

| 数据系列 | 内容 |
|---|---|
| `BF` | FLAG |
| `AP` | 位置（通常存与上一条的差值） |
| `RL` | 读段长度 |
| `MQ` | MAPQ |
| `RN` | 读段名 |
| `QS` | 碱基质量 |
| `FN` / `FC` / `FP` | 每条读段与参考的差异：个数、类型、位置 |
| `BS` / `BA` / `IN` | 替换碱基、未参考比对的碱基、插入序列 |
| `NS` / `NP` / `TS` | mate 的参考、位置、模板长度 |
| `TL` 与标签 block | 可选标签 |

同类数据放在一起时规律更明显，压缩率更高；也意味着读取只需要的字段时，可以跳过其他 block。

### 参考序列依赖

读段碱基按“与参考的差异”保存。例如参考 `ACGTACGT` 上一条完全匹配的读段，只需记录起点和长度；有一个错配时，再记录一个位置和替换碱基。恢复 SEQ 必须使用**压缩时的同一份参考**：

- 每个 slice header 保存它覆盖的参考片段的 **MD5**，解码时据此检查参考是否一致。
- SAM header 的 `@SQ M5:` 也记录每条参考序列的 MD5。
- 文件可以把参考片段**嵌入**到 slice 中（embedded reference），或完全不使用参考压缩（no-ref）；这两种情况读取时不需要外部参考。

### MD / NM 标签

CRAM 通常不存 `MD` 和 `NM` 标签，因为它们可以从参考和差异重新计算。htslib 解码时默认生成它们，OmicsWeave 默认同样生成。

## 关键概念

### 哪些字段需要参考

| 请求的字段 | 是否需要外部参考 |
|---|---|
| 固定列（start、end、flag、mapq、mate、tlen） | 不需要 |
| `cigar`、`name`、`qual` | 不需要 |
| `seq` | 需要（嵌入参考或无参考压缩的文件除外） |
| `aux` 中生成的 `MD` / `NM` | 需要 |

在没有参考的情况下请求 `seq`，会报 `InvalidError` 说明需要参考。

### 参考必须一致

不同来源的 hg38（例如 UCSC 与 GATK bundle）在染色体命名、是否含 decoy / alt contig、甚至个别碱基上可能不同。用错参考时 MD5 校验会失败；关闭校验（`verify_checksums=False`）只会跳过检查，得到的碱基可能是错的。

### 多参考 slice

未定位读段或零散比对较多时，一个 slice 可能包含多条参考序列的记录（multi-reference slice）。CRAI 为这样的 slice 中每条参考各写一行。

## OmicsWeave 中的支持

| 功能 | 支持情况 |
|---|---|
| 版本 | CRAM 3.0、3.1；2.x 及更早、4.0 报 unsupported |
| 格式识别 | 按 `CRAM` magic 识别为 `Cram` |
| block 压缩方法 | raw、gzip、bzip2、lzma、rANS 4x8、rANS Nx16、算术编码、fqzcomp、名称分词（tok3） |
| 数据系列编码 | EXTERNAL、HUFFMAN、BETA、SUBEXP、GAMMA、BYTE_ARRAY_LEN、BYTE_ARRAY_STOP |
| 区域查询 | `AlignmentReader.fetch`，记录与 htslib `sam_itr_queryi` 一致，按文件顺序，结果结构与 BAM 相同 |
| 按需解码 | 只解压所请求字段对应的数据系列 |
| 参考 | 普通 FASTA + `.fai`、bgzip FASTA + `.fai` + `.gzi`、2bit，或已打开的 `SequenceReader` |
| 校验 | 默认检查 block CRC32 和 slice 参考 MD5 |
| 建索引 | `Cram.build_index()` 写 CRAI，不需要参考 |
| 缓存与并发 | 可接 `BlockCache`，缓存解析好的 container 和解码后的 slice；可 `clone()`、pickle |
| 转换 | CRAM → BAM，slice 并行解码，可同时建 BAI / CSI |

限制：

- 不写 CRAM，不提供 BAM → CRAM。
- quant 不直接读取 CRAM，先转换为 BAM。
- 不会从网络（如 EBI 参考服务器）自动下载参考，需要提供本地参考文件。

## 准备索引

```python
from omicsweave.formats import Cram

index = Cram("reads.cram").build_index()
print(index)  # Crai('reads.cram.crai')
```

输入需要按坐标排序。建 CRAI 只解码多参考 slice 的定位信息，不需要外部参考，与 `samtools index` 写出的索引相同。详见 [CRAI](crai.md)。

## 区域查询

```python
from omicsweave.query import AlignmentReader

with AlignmentReader("reads.cram", reference="genome.fa.gz") as reads:
    hits = reads.fetch(("chr1", 100, 200), fields="all")
    print(hits.start, hits.end, hits.flag)
    print(hits.seqs.to_list())
```

返回的 `Alignments` 与 BAM 完全相同，字段说明见 [BAM](bam.md)。

### 不带参考读取

```python
with AlignmentReader("reads.cram") as reads:
    hits = reads.fetch(("chr1", 100, 200), fields=("name", "cigar", "qual"))
```

只需要位置、比对结构或质量时可以不提供参考，读取也更快：碱基和差异相关的 block 不会被解压。嵌入参考或无参考压缩的 CRAM 不提供参考也能恢复 `seq`。

### 共享参考 reader

```python
from omicsweave.cache import BlockCache
from omicsweave.query import AlignmentReader, SequenceReader

cache = BlockCache(max_bytes=512 * 1024**2)
genome = SequenceReader("genome.fa.gz", cache=cache)

sample_a = AlignmentReader("a.cram", reference=genome, cache=cache)
sample_b = AlignmentReader("b.cram", reference=genome, cache=cache)
```

传入已打开的 `SequenceReader` 时，多个 CRAM reader 共享它，由调用方负责关闭；传入路径时，reader 自己打开参考（使用同一个 `cache`）。

### 解码选项

| 参数 | 默认 | 作用 |
|---|---|---|
| `reference` | `None` | 压缩时使用的参考 |
| `index` | 自动查找 | 默认依次找 `<cram>.crai`、`<去扩展名>.crai` |
| `verify_checksums` | `True` | 检查 block CRC32 和 slice 参考 MD5 |
| `decode_md` | `True` | 为没有存储 MD / NM 的记录生成这两个标签 |
| `name_prefix` | 文件名 | 没有存储读段名的记录命名为 `<name_prefix>:<n>` |

这些解码选项只适用于 CRAM；对 BAM 传入 `reference` 或这些选项会报错。

## 转换为 BAM

```sh
ow convert reads.cram reads.bam --reference genome.fa.gz --index -t 8
```

```python
from omicsweave import convert

report = convert.run("reads.cram", "reads.bam", reference="genome.fa.gz", index=True)
print(report.records, report.companions)
```

输出 BAM 的数据流与 htslib 读取 CRAM 再写 BAM 的结果一致（`samtools view -b`，不加 `@PG`）：header 由 CRAM header 重建，记录中包含生成的 MD / NM。转换保留记录顺序，因此建索引仍要求 CRAM 已按坐标排序。`--no-verify` 跳过校验。

## 常见问题

| 现象 | 原因与处理 |
|---|---|
| 请求 `seq` 时报需要参考 | 提供压缩时的参考；或只请求不依赖参考的字段 |
| MD5 校验失败 | 参考与压缩时不同；核对参考版本、染色体命名与 `@SQ M5` |
| 缺少 CRAI | 调用 `Cram(path).build_index()` |
| 想用 BAI / CSI | CRAM 只使用 CRAI |
| quant 不接受 CRAM | 先 `ow convert reads.cram reads.bam --reference genome.fa.gz` |
| 报版本不支持 | 仅支持 CRAM 3.0 / 3.1；旧版本可用 samtools 重新写出 |

## 相关格式

- [BAM](bam.md)：查询结果结构与字段说明，quant 的输入格式。
- [CRAI](crai.md)：CRAM 的区域索引。
- [FASTA](fasta.md)、[2bit](2bit.md)：CRAM 解码使用的参考。
- [SAM](sam.md)：比对记录各字段的含义。

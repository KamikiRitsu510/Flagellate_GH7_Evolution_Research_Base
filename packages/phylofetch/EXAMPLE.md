# phylofetch — 实战示例

> **场景**：你在研究白蚁肠道鞭毛虫共生体的 GH7 纤维素酶基因，从文献和 DIAMOND 比对中筛出了 25 个蛋白 accession，需要：① 下载所有 CDS 核苷酸序列用于密码子比对和选择压力分析；② 下载完整分类谱系做元数据表格。全部流程在一个 Python 脚本里完成。

---

## 查看命令帮助

```bash
# 主帮助
phylofetch --help

# 子命令帮助
phylofetch cds --help
phylofetch taxonomy --help
```

`phylofetch --help` 输出示例：

```
usage: phylofetch [-h] [--email EMAIL] [--verbose] [--quiet] [--retries N]
                  COMMAND ...

Fetch CDS sequences and taxonomy lineages from NCBI.

commands:
  cds        Fetch CDS nucleotide sequences.
  taxonomy   Fetch taxonomy lineages.

options:
  -h, --help          show this help message and exit
  --email, -e EMAIL   Valid email address for NCBI Entrez (required).
  --verbose, -v       Print per-accession progress (default: on).
  --quiet, -q         Suppress progress messages.
  --retries, -r N     Retry attempts on network failure (default: 3).

Examples
--------
  phylofetch cds --input accs.txt --output sequences.fasta --email you@example.com
  phylofetch taxonomy --input accs.txt --excel taxonomy.xlsx --email you@example.com
  phylofetch cds BAB64565.3 AAY83390.3 --email you@example.com --output out.fasta
```

---

## 0. 准备 accession 列表

创建 `gh7_accessions.txt`，每行一个：

```
# GH7 cellulase accessions — flagellate symbionts + bacterial refs
# 格式：每行一个 NCBI accession，# 开头为注释

# 鞭毛虫共生体（Pseudotrichonympha grassii）
BAB64565.3
BAB64566.1
BAB64553.1
BAB64554.1
BAB64555.1
BAB64556.1
BAB64557.1
BAB64558.1
BAB64559.1
BAB64560.1
BAB64561.1
BAB64562.1
BAB64563.2
BAB64564.2

# 鞭毛虫（Holomastigotoides mirabile）
AAY83390.3
AAY83391.1
BAC07551.1
BAC07552.1

# 细菌参考序列（HGT 供体候选）
WP_165634028.1
WP_004955475.1
WP_009974570.1

# 真菌外类群
CAD48174.1
XP_006965674.1
```

---

## 1. 命令行快速使用

### 下载 CDS

```bash
phylofetch cds \
  --input gh7_accessions.txt \
  --output GH7_CDS_raw.fasta \
  --email you@institution.edu

# 终端输出：
# [1/23] Fetching BAB64565.3 ... OK (1035 bp)
# [2/23] Fetching BAB64566.1 ... OK (1038 bp)
# [3/23] Fetching BAB64553.1 ... OK (1032 bp)
# ...
# [21/23] Fetching WP_165634028.1 ... OK (975 bp)
# [22/23] Fetching CAD48174.1 ... OK (1122 bp)
# [23/23] Fetching XP_006965674.1 ... OK (1086 bp)
#
# Done: 23/23 sequences retrieved.
# Saved 23 sequences → GH7_CDS_raw.fasta
```

### 下载分类谱系

```bash
phylofetch taxonomy \
  --input gh7_accessions.txt \
  --excel GH7_taxonomy.xlsx \
  --tsv GH7_taxonomy.tsv \
  --email you@institution.edu

# 终端输出：
# [1/23] BAB64565.3 ... Metamonada
# [2/23] BAB64566.1 ... Metamonada
# ...
# Saved → GH7_taxonomy.xlsx
# Saved → GH7_taxonomy.tsv
```

---

## 2. Python 脚本完整流程

```python
# gh7_fetch.py
from phylofetch import set_email, batch_fetch_cds, batch_fetch_taxonomy

# ── 配置 ──────────────────────────────────────────────────────────────────
set_email("you@institution.edu")

ACCESSIONS = [
    # Pseudotrichonympha grassii
    "BAB64565.3", "BAB64566.1", "BAB64553.1", "BAB64554.1",
    "BAB64555.1", "BAB64556.1", "BAB64557.1", "BAB64558.1",
    "BAB64559.1", "BAB64560.1", "BAB64561.1", "BAB64562.1",
    "BAB64563.2", "BAB64564.2",
    # Holomastigotoides mirabile
    "AAY83390.3", "AAY83391.1", "BAC07551.1", "BAC07552.1",
    # 细菌参考
    "WP_165634028.1", "WP_004955475.1", "WP_009974570.1",
    # 外类群
    "CAD48174.1", "XP_006965674.1",
]

# ── 1. 下载 CDS ───────────────────────────────────────────────────────────
print("=== 下载 CDS 序列 ===")
seqs = batch_fetch_cds(
    ACCESSIONS,
    output_fasta="GH7_CDS_raw.fasta",
    retries=3,
    verbose=True,
)

print(f"\n成功：{len(seqs)}/{len(ACCESSIONS)} 条")
print(f"序列长度范围：{min(len(s) for s in seqs.values())}–"
      f"{max(len(s) for s in seqs.values())} bp")

# 找出下载失败的
failed = [a for a in ACCESSIONS if a not in seqs]
if failed:
    print(f"⚠ 未找到：{failed}")

# ── 2. 下载分类谱系 ───────────────────────────────────────────────────────
print("\n=== 下载分类信息 ===")
df = batch_fetch_taxonomy(
    ACCESSIONS,
    output_excel="GH7_taxonomy.xlsx",
    output_tsv="GH7_taxonomy.tsv",
    extra_ranks=["subphylum"],   # 额外获取亚门级别
    verbose=True,
)

print("\n=== 分类摘要 ===")
print(df[["accession", "superkingdom", "phylum", "class", "genus"]].to_string())
```

运行：

```bash
python gh7_fetch.py
```

---

## 3. 查看结果

### GH7_CDS_raw.fasta 内容（前几行）

```
>BAB64565.3
ATGAAAAATTTGCTTATTTTGTTCTTGATCATAGCAGTCACAGCAAGCAGAAATCCTGCAGATCCAAGT
CCTGATACCAAAAGCAAAGATAATGACAGCAGCCCGTACATATAAGCCTGATATTACAGAAGACATTATT
...
>BAB64566.1
ATGAAAAATTTGCTTATTTTCTTGATCATAGCAATCACAGCAAGCAGAAATCCTGCAGATCCAAGTCCT
...
```

### GH7_taxonomy.tsv 内容（前几行）

```
accession       taxid   superkingdom  phylum      class           genus
BAB64565.3      5722    Eukaryota     Metamonada  Parabasalia     Pseudotrichonympha
BAB64566.1      5722    Eukaryota     Metamonada  Parabasalia     Pseudotrichonympha
BAB64553.1      5722    Eukaryota     Metamonada  Parabasalia     Pseudotrichonympha
AAY83390.3      286380  Eukaryota     Metamonada  Parabasalia     Holomastigotoides
WP_165634028.1  1505    Bacteria      Firmicutes  Clostridia      Pseudobacteroides
CAD48174.1      34608   Eukaryota     Chytridiomycota  Neocallimastigomycetes  Piromyces
```

---

## 4. 后续对接 msaclean（典型工作流）

```bash
# 1. 用 MAFFT 比对下载的 CDS（或用 pal2nal 做密码子比对）
mafft --localpair --maxiterate 1000 GH7_CDS_raw.fasta > GH7_CDS_aligned.fasta

# 2. 用 msaclean 清洗比对，保护核心鞭毛虫序列
msaclean clean GH7_CDS_aligned.fasta \
  -o GH7_CDS_clean.fasta \
  -w core_accessions.txt \
  --max-gap 0.9
```

---

## 注意事项

- **NCBI 限速**：默认 0.4 秒/请求（≈2.5 req/s），23 条序列约需 20–30 秒。有 API key 可在 `set_email()` 后设 `from Bio import Entrez; Entrez.api_key = "your_key"` 加速。
- **蛋白 accession 取 CDS**：`BAB64565.3` 是蛋白 accession，程序会从 GenBank 蛋白记录的 `coded_by` 字段提取核酸坐标再去取序列，多一次网络请求，正常现象。
- **序列找不到**：部分较老的 accession 可能已被 NCBI suppress，函数会打印警告并跳过，不会中断整体任务。

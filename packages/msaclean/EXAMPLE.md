# msaclean — 实战示例

> **场景**：你用 CD-HIT 对 NCBI 检索到的 800 余条 GH7 纤维素酶蛋白序列做了聚类（90% identity），保留代表序列 187 条。用 MAFFT L-INS-i 全局比对后发现比对质量很差——大量序列几乎全是 gap，还有几条远缘序列把整个比对拉得很长。你有 17 条来自鞭毛虫共生体的核心序列必须保留。目标：清洗出一个高质量的 ~50 条序列比对，用于后续 IQ-TREE 建树。

---

## 查看命令帮助

```bash
# 主帮助
msaclean --help

# 子命令帮助
msaclean clean --help
msaclean iterative --help
msaclean stats --help
```

`msaclean --help` 输出示例：

```
usage: msaclean [-h] [--verbose] [--quiet] COMMAND ...

Clean and curate multiple sequence alignments.

commands:
  clean       Remove empty/gappy sequences from an existing alignment.
  iterative   Iteratively realign (MAFFT) and remove gappiest sequences.
  stats       Report gap statistics for each sequence in an alignment.

options:
  -h, --help   show this help message and exit
  --verbose    Verbose output.
  --quiet      Suppress progress messages.

Examples
--------
  msaclean clean GH7_aligned.fasta -o GH7_clean.fasta -w whitelist.txt
  msaclean stats GH7_aligned.fasta --top 20 -w whitelist.txt
  msaclean iterative GH7_raw.fasta -o GH7_curated.fasta \
      -w whitelist.txt --stop-at 20 --rounds 2
```

`msaclean iterative --help` 输出示例：

```
usage: msaclean iterative [-h] [--output FILE] [--whitelist FILE]
                           [--stop-at N] [--rounds N] [--max-rounds N]
                           [--max-gap FRAC] [--mafft-threads N]
                           [--mafft-mode {localpair,globalpair,genafpair,auto}]
                           [--no-trimal]
                           input

positional arguments:
  input                 Input FASTA (unaligned or pre-aligned).

options:
  -h, --help            show this help message and exit
  --output, -o FILE     Output alignment FASTA.
  --whitelist, -w FILE  Text file with accessions to never remove.
  --stop-at N           Stop when this many sequences remain.
  --rounds N            Sequences to remove per round (default 2).
  --max-rounds N        Maximum number of rounds (default 20).
  --max-gap FRAC        Only remove sequences above this gap fraction.
  --mafft-threads N     Threads for MAFFT (default 4).
  --mafft-mode ...      MAFFT alignment mode (default: localpair).
  --no-trimal           Skip trimAl column-trimming step.
```

---

## 0. 准备白名单文件

创建 `whitelist.txt`（这些序列无论 gap 多少都不会被删除）：

```
# GH7 核心鞭毛虫序列白名单
# Pseudotrichonympha grassii
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
# Holomastigotoides mirabile
AAY83390.3
AAY83391.1
BAC07551.1
```

---

## 第一步：先看看比对有多糟糕

```bash
msaclean stats GH7_aligned_187.fasta --top 25 -w whitelist.txt
```

输出：

```
Alignment: 187 sequences × 3842 columns

Accession                       Gap%   Residues   Protected
--------------------------------------------------------------
WP_098765432.1                 98.2%         69
OXZ87654.1                     97.6%         92
KAF23456.1                     96.9%        119
RZX12345.1                     95.4%        176
...
BAB64565.3                     43.2%       2183           ✓
BAB64553.1                     41.8%       2237           ✓
AAY83390.3                     44.6%       2130           ✓
...
```

结论：约 30 条序列 gap > 90%，需要清除。白名单序列 gap 在 40–50% 属于正常（比对柱很长，有合理 gap）。

---

## 第二步：一键清洗全空序列

```bash
msaclean clean GH7_aligned_187.fasta \
  -o GH7_aligned_clean.fasta \
  -w whitelist.txt \
  --max-gap 0.90
```

输出：

```
=== msaclean report ===
Input sequences : 187
Protected       : 17
Removed (empty) : 0
Removed (gappy) : 31
Kept            : 156
Output          : GH7_aligned_clean.fasta

Gappy sequences : WP_098765432.1, OXZ87654.1, KAF23456.1, ...
```

156 条，从 187 降到 156。但比对质量还不够好，需要继续迭代。

---

## 第三步：迭代清洗到目标数量

```bash
msaclean iterative GH7_aligned_clean.fasta \
  -o GH7_curated_final.fasta \
  -w whitelist.txt \
  --stop-at 50 \
  --rounds 3 \
  --mafft-mode localpair \
  --mafft-threads 8
```

运行过程输出：

```
Starting with 156 sequences (17 whitelisted).
Round  0:  156 seqs,   2184 columns  → removing: RZX12345.1, KGH23456.1, MNP98765.1
Round  1:  153 seqs,   2201 columns  → removing: OPQ34567.1, LMN45678.1, QRS56789.1
Round  2:  150 seqs,   2198 columns  → removing: TUV67890.1, WXY78901.1, ZAB89012.1
Round  3:  147 seqs,   2176 columns  → removing: BCD90123.1, EFG01234.1, HIJ12345.1
...
Round 27:   53 seqs,   1843 columns  → removing: KLM23456.1, NOP34567.1, QRS45678.1
Round 28:   50 seqs,   1821 columns  → STOP (target 50 reached)

Round   Seqs  Aln len  Removed
     0   156     2184  RZX12345.1, KGH23456.1, MNP98765.1
     1   153     2201  OPQ34567.1, LMN45678.1, QRS56789.1
  ...
    28    50     1821  KLM23456.1, NOP34567.1, QRS45678.1
Stopped: reached stop_at_n=50

Final: 50 sequences → GH7_curated_final.fasta
```

---

## 第四步：Python 脚本版（更灵活）

```python
# curate_gh7.py
from Bio import AlignIO
from msaclean import clean_alignment_file, gappiest_sequences, iterative_clean

WHITELIST = {
    "BAB64565.3", "BAB64566.1", "BAB64553.1", "BAB64554.1",
    "BAB64555.1", "BAB64556.1", "BAB64557.1", "BAB64558.1",
    "BAB64559.1", "BAB64560.1", "BAB64561.1", "BAB64562.1",
    "BAB64563.2", "BAB64564.2",
    "AAY83390.3", "AAY83391.1", "BAC07551.1",
}

# ── 第一遍：快速清洗全空序列 ──────────────────────────────────────────────
report = clean_alignment_file(
    "GH7_aligned_187.fasta",
    "GH7_step1_clean.fasta",
    whitelist=WHITELIST,
    max_gap_fraction=0.90,
)
print(f"第一遍清洗：{report['original']} → {report['kept']} 条")
print(f"  删除 {report['removed_gappy']} 条 gap>90% 序列")

# ── 第二遍：看看还有哪些序列拖后腿 ──────────────────────────────────────
aln = AlignIO.read("GH7_step1_clean.fasta", "fasta")
worst = gappiest_sequences(aln, n=10, whitelist=WHITELIST)
print("\n最差的 10 条（gap 最多）：")
for acc, gf in worst:
    print(f"  {acc:30s}  {gf:.1%}")

# ── 第三遍：迭代清洗 ──────────────────────────────────────────────────────
result = iterative_clean(
    "GH7_step1_clean.fasta",
    output_fasta="GH7_curated_final.fasta",
    whitelist=WHITELIST,
    remove_per_round=3,
    stop_at_n=50,
    max_rounds=30,
    mafft_threads=8,
    mafft_mode="localpair",
    use_trimal=True,
    verbose=True,
)

print("\n=== 迭代清洗完成 ===")
print(result.log_table)
print(f"\n最终序列数：{len(result.final_alignment)}")
print(f"比对长度：{result.final_alignment.get_alignment_length()} 列")
```

---

## 最终结果验证

```bash
# 用 msaclean stats 验证最终比对质量
msaclean stats GH7_curated_final.fasta -w whitelist.txt

# 期望输出：
# Alignment: 50 sequences × 1821 columns
#
# Accession                       Gap%   Residues   Protected
# ------------------------------------------------------------------
# PXY12345.1                     62.1%        691
# OAB56789.1                     58.4%        757
# BAB64565.3                     41.2%       1071           ✓
# AAY83390.3                     42.8%       1042           ✓
# ...
# （所有序列 gap < 65%，质量显著提升）
```

```bash
# 对接 IQ-TREE 建树
iqtree2 -s GH7_curated_final.fasta \
  -m TEST -bb 1000 -nt AUTO \
  --prefix GH7_final
```

---

## 注意事项

- **白名单匹配**：程序取 FASTA header 第一个空格前的字符串作为 accession。`>BAB64565.3 Pseudotrichonympha grassii` 会匹配到 `BAB64565.3`。
- **迭代轮数估算**：从 N 条序列清洗到 M 条，每轮删 k 条，大概需要 (N-M)/k 轮。加上每轮 MAFFT L-INS-i 的时间，100 条序列约 2–5 分钟/轮。序列太多可先用 `--mafft-mode auto` 快速预筛。
- **trimAl 没装**：加 `--no-trimal` 跳过列修剪，照样能跑，只是每轮比对不做自动列过滤。
- **结果不稳定**：迭代清洗的结果依赖初始顺序，如果对结果有疑虑，可以多跑几次或调整 `--rounds` 参数对比。

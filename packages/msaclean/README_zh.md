# msaclean

**用于系统发育分析的多序列比对清洗工具包。**

本包将 `clean_with_whitelist.py` 和 `iterative_clean.sh` 的核心逻辑提炼为可复用的 Python 包，解决比对后序列质量控制中最常见的两类问题：

- **全空序列**：聚类后的序列比对进来几乎全是 gap，对树没有贡献。
- **迭代清洗**：某些序列严重拉低比对质量，需要一轮轮删掉最差的，每轮重新比对。

白名单（whitelist）机制贯穿全包：指定的 accession 无论 gap 比例多高都不会被删除。

---

## 安装

```bash
pip install msaclean
# 或从本地安装
pip install path/to/packages/msaclean
```

**迭代清洗功能额外需要**（需在系统 PATH 中）：

```bash
conda install -c bioconda mafft trimal
```

---

## 快速上手

### 场景一：清洗已有比对文件（最常用）

```python
from msaclean import clean_alignment_file

WHITELIST = {
    "BAB64565.3", "AAY83390.3", "BAB64553.1", "BAB64554.1",
    "BAB64555.1", "BAB64556.1", "BAB64557.1", "BAB64558.1",
}

report = clean_alignment_file(
    "GH7_aligned.fasta",    # 输入：已比对的 FASTA
    "GH7_clean.fasta",      # 输出
    whitelist=WHITELIST,
    max_gap_fraction=0.95,  # 删除 gap 比例超过 95% 的序列
)

print(f"原始：{report['original']} 条")
print(f"保留：{report['kept']} 条")
print(f"删除（全空）：{report['removed_empty']} 条")
print(f"删除（过多gap）：{report['removed_gappy']} 条")
print(f"白名单保护：{report['protected']} 条")
```

### 场景二：在内存中操作 Biopython 对象

```python
from Bio import AlignIO
from msaclean import clean_alignment

aln = AlignIO.read("GH7_aligned.fasta", "fasta")
clean, report = clean_alignment(aln, whitelist=WHITELIST, max_gap_fraction=0.9)

# 继续处理 clean（MultipleSeqAlignment 对象）
AlignIO.write(clean, "GH7_clean.fasta", "fasta")
```

### 场景三：先查看哪些序列最差

```python
from Bio import AlignIO
from msaclean import gappiest_sequences

aln = AlignIO.read("GH7_aligned.fasta", "fasta")
worst = gappiest_sequences(aln, n=10, whitelist=WHITELIST)

for acc, gap_frac in worst:
    print(f"{acc:30s}  {gap_frac:.1%} gaps")
```

先查看再决定阈值，比直接乱设 `max_gap_fraction` 更稳妥。

### 场景四：迭代清洗（需要 MAFFT）

```python
from msaclean import iterative_clean

result = iterative_clean(
    "GH7_clustered.fasta",       # 输入：未比对或已比对均可
    output_fasta="GH7_curated.fasta",
    whitelist=WHITELIST,
    remove_per_round=2,          # 每轮删除最差的 2 条
    stop_at_n=20,                # 剩 20 条时停止
    mafft_mode="localpair",      # L-INS-i 模式，精度最高
    use_trimal=True,             # 每轮比对后用 trimAl 修剪列
)

print(result.log_table)
# Round   Seqs  Aln len  Removed
#     0    187     1254  XP_034521.1, WP_165634028.1
#     1    185     1189  ...
# Stopped: reached stop_at_n=20
```

---

## 命令行用法

```bash
# 查看比对文件的 gap 统计（不修改文件）
msaclean stats GH7_aligned.fasta --top 20 -w whitelist.txt

# 清洗：删除全空和 gap>90% 的序列
msaclean clean GH7_aligned.fasta -o GH7_clean.fasta \
    -w whitelist.txt --max-gap 0.9

# 迭代清洗：每轮删 2 条，剩 20 条停止
msaclean iterative GH7_raw.fasta -o GH7_curated.fasta \
    -w whitelist.txt --stop-at 20 --rounds 2

# 用更快的 MAFFT 模式，不用 trimAl
msaclean iterative GH7_raw.fasta -o GH7_curated.fasta \
    --mafft-mode auto --no-trimal --stop-at 30
```

`whitelist.txt` 格式：每行一个 accession，`#` 开头为注释。

---

## 注意事项

**白名单 accession 的匹配逻辑**：本包从序列 ID 中取**第一个空格前**的字符串作为 accession（即 `record.id.split()[0]`）。FASTA header 格式为 `>BAB64565.3 Pseudotrichonympha ...` 时会正确提取 `BAB64565.3`。如果你的 header 格式不同，请先用 Biopython 手动修改 `record.id`。

**迭代清洗的输入**：`iterative_clean()` 接受未比对或已比对的 FASTA。每轮开始前都会用 MAFFT 重新比对，所以输入是否已比对不影响结果。

**trimAl 不可用时**：传入 `use_trimal=False` 跳过列修剪步骤；或命令行加 `--no-trimal`。没有 trimAl 时清洗照样能跑，只是每轮比对不做列过滤。

**内存占用**：`iterative_clean()` 会把每轮的中间结果写到系统临时目录（`/tmp/msaclean_*/`），任务结束后自动清理。序列数超过 500 条时 L-INS-i 会很慢，建议先用 `--mafft-mode auto`。

**只想删空序列不想动其他的**：设 `max_gap_fraction=1.0`（默认值）即可，只删除残基数为 0 的全空序列，其他不动。

---

## 函数速查

| 函数/类 | 说明 |
|---------|------|
| `clean_alignment_file(input, output, ...)` | 文件级清洗，最常用 |
| `clean_alignment(aln, ...)` | 对象级清洗，返回 `(alignment, report)` |
| `filter_empty_seqs(aln, whitelist, min_residues)` | 只删全空/近空序列 |
| `filter_gappy_seqs(aln, whitelist, max_gap_fraction)` | 只删超阈值 gap 序列 |
| `gappiest_sequences(aln, n, whitelist)` | 返回 gap 最多的 n 条 accession |
| `protect_sequences(aln, whitelist)` | 将比对拆分为保护/非保护两部分 |
| `iterative_clean(input_fasta, ...)` | 迭代重比对+清洗（需 MAFFT） |
| `IterativeCleanResult.log_table` | 输出轮次统计表格的字符串 |

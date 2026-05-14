# phylofetch

**从 NCBI 批量下载 CDS 序列和分类学谱系信息的 Python 工具包。**

本包把 Biopython 的 Entrez 接口封装成限速、自动重试、批量友好的函数，专为系统发育研究中"给一堆 accession 号批量取序列/取分类"的场景而生。原始脚本来自 GH7 纤维素酶 HGT 研究，现提炼为通用包供复用。

---

## 安装

```bash
pip install phylofetch
# 或从本地安装
pip install path/to/packages/phylofetch
```

可选依赖（ete3 树操作）：

```bash
pip install "phylofetch[full]"
```

---

## 使用前必读：NCBI 邮箱设置

NCBI 的 Entrez 接口**要求**每次请求都附带一个合法邮箱地址，否则会被拒绝服务。请在脚本最开头调用：

```python
from phylofetch import set_email
set_email("you@example.com")   # 换成你自己的邮箱
```

> **注意**：不需要注册，填任何真实邮箱即可。漏掉这一步会直接抛出 `RuntimeError`。

---

## 核心功能

### 1. 获取单条 CDS 序列

```python
from phylofetch import set_email, fetch_cds

set_email("you@example.com")
seq = fetch_cds("BAB64565.3")
print(f"序列长度：{len(seq)} bp")
print(seq[:60])
```

`fetch_cds()` 内部逻辑：
1. 先在 **核酸数据库** 搜索该 accession
2. 找不到就去 **蛋白数据库**，从 `coded_by` qualifier 解析坐标再取核酸子序列
3. 返回纯字符串（ACGT），找不到时返回 `None`

### 2. 批量获取 CDS，保存为 FASTA

```python
from phylofetch import batch_fetch_cds

accs = ["BAB64565.3", "AAY83390.3", "BAB64553.1", "XP_006965674.1"]
seqs = batch_fetch_cds(accs, output_fasta="my_CDS.fasta", verbose=True)

print(f"成功下载：{len(seqs)}/{len(accs)} 条")
# seqs 是字典：{accession: 序列字符串}
```

运行时会逐条打印进度，例如：

```
[1/4] Fetching BAB64565.3 ... OK (1035 bp)
[2/4] Fetching AAY83390.3 ... OK (1038 bp)
...
Saved 4 sequences → my_CDS.fasta
```

### 3. 获取单条分类信息

```python
from phylofetch import fetch_taxonomy

tax = fetch_taxonomy("BAB64565.3")
print(tax["phylum"], "/", tax["class"])
# Metamonada / Parabasalia

# 返回的字典包含：
# accession, taxid, lineage_string,
# superkingdom, kingdom, phylum, class, order, family, genus, species
# 缺少某个级别时值为 None
```

### 4. 批量获取分类，保存为 Excel / TSV

```python
from phylofetch import batch_fetch_taxonomy

accs = ["BAB64565.3", "AAY83390.3", "WP_165634028.1"]
df = batch_fetch_taxonomy(
    accs,
    output_excel="taxonomy.xlsx",
    output_tsv="taxonomy.tsv",
    extra_ranks=["subphylum"],   # 额外获取的分类级别（可选）
)
print(df[["accession", "phylum", "class", "genus"]])
```

---

## 命令行用法

安装后直接在终端使用：

```bash
# 从文件读取 accession 列表，下载 CDS
phylofetch cds --input accs.txt --output sequences.fasta --email you@example.com

# 直接在命令行写 accession
phylofetch cds BAB64565.3 AAY83390.3 --email you@example.com

# 获取分类，同时保存 Excel 和 TSV
phylofetch taxonomy --input accs.txt \
    --excel taxonomy.xlsx --tsv taxonomy.tsv \
    --email you@example.com

# 包含额外分类级别
phylofetch taxonomy --input accs.txt \
    --extra-ranks subphylum tribe \
    --tsv full_taxonomy.tsv --email you@example.com
```

`accs.txt` 格式：每行一个 accession，`#` 开头的行视为注释忽略。

---

## 注意事项

**限速**：NCBI 对未注册用户限制每秒 ≤3 次请求。本包内置 0.4 秒最小间隔，无需手动控制。如有 API key（`Entrez.api_key = "your_key"`），速率可提升至 10 次/秒，但包本身不会自动加速，需自行调整 `RateLimiter`。

**重试**：所有请求在网络错误时自动重试（默认 3 次，指数退避），适合不稳定的网络环境。`retries` 参数可调整次数。

**蛋白 vs 核酸 accession**：两种都支持。蛋白 accession（如 `BAB64565.3`）会自动从 `coded_by` 字段提取 CDS 坐标；核酸 accession 直接取全长序列。

**序列找不到**：`fetch_cds()` 返回 `None`，`batch_fetch_cds()` 跳过该条目（不会中断批量任务）；错误信息会写入日志。

---

## 函数速查

| 函数 | 说明 |
|------|------|
| `set_email(email)` | 配置 NCBI 邮箱（必须最先调用） |
| `fetch_cds(accession)` | 获取单条 CDS，返回字符串或 None |
| `batch_fetch_cds(accs, output_fasta=...)` | 批量获取 CDS，返回字典 |
| `fetch_taxonomy(accession)` | 获取单条分类，返回字典 |
| `batch_fetch_taxonomy(accs, output_excel=..., output_tsv=...)` | 批量获取分类，返回 DataFrame |
| `RateLimiter(min_interval)` | 手动创建限速器（高级用法） |

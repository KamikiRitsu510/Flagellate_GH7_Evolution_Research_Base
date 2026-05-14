#!/bin/bash
echo "======================================"
echo "  从 NCBI 获取物种名称（改进版）"
echo "======================================"

EMAIL="kamikiritsu510@gmail.com"
FASTA_FILE="data/final_matrix.fasta"
OUTPUT_FILE="accession_species.txt"

# 检查依赖
if ! command -v python3 &> /dev/null; then
    echo "❌ Python3 未安装"
    exit 1
fi
python3 -c "import Bio" 2>/dev/null || pip3 install biopython --quiet

# 提取登录号
grep "^>" "$FASTA_FILE" | sed 's/^>//' | awk '{print $1}' > acc_list.txt
TOTAL=$(wc -l < acc_list.txt | xargs)
echo "✅ 共 $TOTAL 个登录号"

# Python 脚本：从 GenBank DEFINITION 提取物种名
cat > fetch_species.py << 'PYEOF'
import sys, time
from Bio import Entrez

Entrez.email = sys.argv[1]

with open("acc_list.txt") as f:
    accessions = [line.strip() for line in f if line.strip()]

results = []
for idx, acc in enumerate(accessions, 1):
    print(f"  [{idx}/{len(accessions)}] {acc} ...", end=" ", flush=True)
    try:
        # 获取 GenBank 记录
        handle = Entrez.efetch(db="protein", id=acc, rettype="gb", retmode="text")
        record = handle.read()
        handle.close()
        # 提取 DEFINITION 行（通常在开头）
        for line in record.split("\n"):
            if line.startswith("DEFINITION"):
                # 示例: "DEFINITION  hypothetical protein [uncultured bacterium]"
                # 物种名出现在方括号内或者最后的部分
                import re
                # 尝试提取 [物种名]
                match = re.search(r'\[([^\]]+)\]', line)
                if match:
                    species = match.group(1)
                else:
                    # 没有方括号，取 DEFINITION 后半部分（通常“ 物种名”）
                    parts = line.split("  ", 1)
                    if len(parts) > 1:
                        species = parts[1].split(",")[0].strip()
                    else:
                        species = "Unknown"
                print(f"→ {species}")
                results.append((acc, species))
                break
        else:
            print("→ Not found")
            results.append((acc, "Not found"))
        time.sleep(0.35)  # 避免请求过快
    except Exception as e:
        print(f"→ Error: {e}")
        results.append((acc, "Error"))
        time.sleep(1)

with open("accession_species.txt", "w") as f:
    f.write("accession\tspecies\n")
    for acc, sp in results:
        f.write(f"{acc}\t{sp}\n")
print("\n✅ 完成！")
PYEOF

python3 fetch_species.py "$EMAIL"

rm -f acc_list.txt fetch_species.py
echo "📄 结果保存到: $OUTPUT_FILE"
head -10 "$OUTPUT_FILE"

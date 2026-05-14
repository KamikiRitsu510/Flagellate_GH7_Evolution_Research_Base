from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

WHITELIST = {
    "AAY83390.3", "AAY83391.1",
    "BAB64553.1", "BAB64554.1", "BAB64555.1", "BAB64556.1",
    "BAB64557.1", "BAB64558.1", "BAB64559.1", "BAB64560.1",
    "BAB64561.1", "BAB64562.1", "BAB64563.2", "BAB64564.2",
    "BAB64565.3", "BAC07551.1", "BAC07552.1"
}

original_alignment = AlignIO.read("GH7_aligned_linsi_final.fasta", "fasta")

clean_records = []
protected_count = 0
removed_empty = 0
skipped_bad = 0

for record in original_alignment:
    try:
        acc = record.id.strip().split()[0]
    except (AttributeError, IndexError):
        # 头部信息损坏，直接跳过
        skipped_bad += 1
        print(f"  跳过损坏记录: {record.id}")
        continue

    # 白名单保护
    if acc in WHITELIST:
        clean_records.append(record)
        protected_count += 1
        continue

    # 非白名单：检查是否全是gap
    residues = str(record.seq).replace('-', '').replace('.', '')
    if len(residues) > 0:
        clean_records.append(record)
    else:
        removed_empty += 1
        print(f"  删除空序列: {record.id}")

clean_alignment = MultipleSeqAlignment(clean_records)
AlignIO.write(clean_alignment, "GH7_linsi_whitelisted.fasta", "fasta")

print(f"\n=== 清洗报告 ===")
print(f"原始序列总数: {len(original_alignment)}")
print(f"白名单保护: {protected_count}/{len(WHITELIST)}")
print(f"删除空序列: {removed_empty}")
print(f"跳过损坏记录: {skipped_bad}")
print(f"保留序列总数: {len(clean_records)}")
print(f"输出文件: GH7_linsi_whitelisted.fasta")

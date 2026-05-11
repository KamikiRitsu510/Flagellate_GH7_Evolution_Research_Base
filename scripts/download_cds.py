#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

Entrez.email = "kamikiritsu@gmail.com"

def get_cds_from_genbank(accession):
    try:
        handle = Entrez.esearch(db="nucleotide", term=accession)
        record = Entrez.read(handle)
        handle.close()
        if record["Count"] == "0":
            handle = Entrez.esearch(db="protein", term=accession)
            rec = Entrez.read(handle)
            handle.close()
            if rec["Count"] == "0":
                return None
            prot_handle = Entrez.efetch(db="protein", id=rec["IdList"][0], rettype="gb", retmode="text")
            prot_record = SeqIO.read(prot_handle, "genbank")
            prot_handle.close()
            for feature in prot_record.features:
                if feature.type == "CDS":
                    if "coded_by" in feature.qualifiers:
                        coded_by = feature.qualifiers["coded_by"][0]
                        match = re.search(r"([A-Z0-9_.]+):(\d+)-(\d+)", coded_by)
                        if match:
                            nuc_acc = match.group(1)
                            start = int(match.group(2))
                            end = int(match.group(3))
                            nuc_handle = Entrez.efetch(db="nucleotide", id=nuc_acc, rettype="fasta", retmode="text", seq_start=start, seq_stop=end)
                            seq_record = SeqIO.read(nuc_handle, "fasta")
                            nuc_handle.close()
                            return str(seq_record.seq)
            return None
        else:
            nuc_handle = Entrez.efetch(db="nucleotide", id=record["IdList"][0], rettype="fasta", retmode="text")
            seq_record = SeqIO.read(nuc_handle, "fasta")
            nuc_handle.close()
            return str(seq_record.seq)
    except Exception as e:
        print(f"  Error for {accession}: {e}")
        return None

def main():
    target_accs = [
        "AAY83390.3", "AAY83391.1", "BAB64553.1", "BAB64565.3",
        "BAC07551.1", "BAC07552.1", "XP_006965674.1",
        "XP_003054277.1", "AGT15838.1"
    ]
    tree_file = "results/GH7_final_534.contree"
    try:
        with open(tree_file) as f:
            content = f.read()
        tips = re.findall(r'([A-Za-z0-9_.]+)(?=[:,)]|$)', content)
        daphnia_tips = [t for t in tips if 'daphnia' in t.lower() or 'pulex' in t.lower()]
        if daphnia_tips:
            print(f"发现水蚤序列: {daphnia_tips[:3]}")
            target_accs.extend(daphnia_tips[:2])
        else:
            print("未发现水蚤序列，将只使用现有登录号。")
    except FileNotFoundError:
        print(f"警告: 树文件 {tree_file} 不存在，跳过添加水蚤序列。")
    
    target_accs = list(dict.fromkeys(target_accs))
    print(f"总共需要下载 {len(target_accs)} 条序列的 CDS:")
    for acc in target_accs:
        print(f"  {acc}")
    
    records = []
    for acc in target_accs:
        print(f"正在下载 {acc} ...")
        cds_seq = get_cds_from_genbank(acc)
        if cds_seq:
            records.append(SeqRecord(Seq(cds_seq), id=acc, description=""))
            print(f"   成功，长度 {len(cds_seq)} bp")
        else:
            print(f"   失败，未找到 CDS")
    
    if records:
        output_file = "raw_cds_subset.fasta"
        with open(output_file, "w") as out:
            SeqIO.write(records, out, "fasta")
        print(f"\n已保存 {len(records)} 条 CDS 序列到 {output_file}")
    else:
        print("没有成功下载任何序列，请检查网络或登录号。")

if __name__ == "__main__":
    main()

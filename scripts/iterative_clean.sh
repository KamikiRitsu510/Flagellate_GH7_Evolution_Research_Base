#!/bin/bash
WORK_DIR="/Users/kamikiritsu/Desktop/unsorted/GH7 origin"
INPUT_FASTA="$WORK_DIR/02_GH7_cdhit_clustered_clean.fasta"
LOG_FILE="$WORK_DIR/iterative_cleaning_log.txt"
CORE_SEQS="AAY83390.3|AAY83391.1|BAB64553.1|BAB64554.1|BAB64555.1|BAB64556.1|BAB64557.1|BAB64558.1|BAB64559.1|BAB64560.1|BAB64561.1|BAB64562.1|BAB64563.2|BAB64564.2|BAB64565.3|BAC07551.1|BAC07552.1"

cd "$WORK_DIR"
cp "$INPUT_FASTA" current_set.fasta
echo "Round,SeqCount,AlignedLength" > "$LOG_FILE"

for round in {0..12}; do
    mafft --localpair --maxiterate 1000 --thread 4 current_set.fasta > current_aln.fasta 2>/dev/null
    trimal -in current_aln.fasta -out current_trimmed.fasta -automated1 2>/dev/null
    
    seq_count=$(grep -c ">" current_trimmed.fasta)
    # 用测序长度近似替代
    aligned_len=$(head -2 current_trimmed.fasta | tail -1 | awk '{print length}')
    
    echo "$round,$seq_count,$aligned_len" >> "$LOG_FILE"
    echo "Round $round: $seq_count sequences, ~$aligned_len columns"
    
    if [ "$seq_count" -le 20 ]; then
        echo "已接近核心序列数量，停止迭代。"
        break
    fi
    
    # 找出gap最多的2条（排除核心序列）
    grep -vE "$CORE_SEQS" current_trimmed.fasta | \
        awk '/^>/ {name=$0; next} {print gsub("-",""), name}' | \
        sort -nr | head -2 | awk '{print $2}' | sed 's/>//' > worst_2.txt
    
    to_remove=$(cat worst_2.txt)
    if [ -z "$to_remove" ]; then
        echo "没有可删除的非核心序列，停止。"
        break
    fi
    
    for id in $to_remove; do
        echo "删除: $id"
        seqkit grep -v -p "$id" current_set.fasta > tmp.fasta
        mv tmp.fasta current_set.fasta
    done
done

echo "迭代清洗完成。日志: $LOG_FILE"

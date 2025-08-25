#!/bin/bash

# 检查输入参数
if [ "$#" -lt 3 ] || [ "$#" -gt 5 ]; then
    echo "Usage: $0 <hmm_directory> <fasta_directory> <output_directory> [hmm_threads=32] [parallel_jobs=24]"
    exit 1
fi

# 获取输入参数
HMM_DIR="$1"
FASTA_DIR="$2"
OUTPUT_DIR="$3"
HMM_THREADS="${4:-8}"       # hmmsearch --cpu 默认为 32
PARALLEL_JOBS="${5:-32}"     # parallel -j 默认为 24

# 检查目录存在性
if [ ! -d "$HMM_DIR" ]; then
    echo "❌ Error: HMM directory '$HMM_DIR' not found."
    exit 1
fi

if [ ! -d "$FASTA_DIR" ]; then
    echo "❌ Error: FASTA directory '$FASTA_DIR' not found."
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 定义并导出函数
process_hmm_fasta() {
    HMM_FILE="$1"
    FASTA_FILE="$2"
    OUTPUT_DIR="$3"
    HMM_THREADS="$4"

    HMMNAME=$(basename "$HMM_FILE" .hmm)
    BASENAME=$(basename "$FASTA_FILE")
    BASENAME=${BASENAME%.fasta}
    BASENAME=${BASENAME%.faa}
    BASENAME=${BASENAME%.gz}

    OUTPUT_PATH="$OUTPUT_DIR/$HMMNAME"
    mkdir -p "$OUTPUT_PATH"
    LOG_FILE="$OUTPUT_PATH/hmmsearch_errors.log"

    tblout_file="$OUTPUT_PATH/${BASENAME}.tblout"
    domtblout_file="$OUTPUT_PATH/${BASENAME}.domtblout"
    temp_file="$OUTPUT_PATH/${BASENAME}_cleaned.faa"

    # 如果输出文件已存在且不为空，跳过
    if [ -s "$tblout_file" ] && [ -s "$domtblout_file" ]; then
        echo "✅ Skipping $domtblout_file (already exists)"
        return
    fi

    # 解压和清理序列
    if [[ "$FASTA_FILE" == *.gz ]]; then
        echo "🧼 Decompressing and cleaning $FASTA_FILE"
        zcat "$FASTA_FILE" | sed '/^>/! s/-//g' > "$temp_file"
    else
        echo "🧼 Cleaning $FASTA_FILE"
        sed '/^>/! s/-//g' "$FASTA_FILE" > "$temp_file"
    fi

    # 执行 hmmsearch
    echo "🚀 Running hmmsearch on $FASTA_FILE with $HMM_FILE"
    /usr/bin/hmmsearch --cpu "$HMM_THREADS" --domtblout "$domtblout_file" --tblout "$tblout_file" "$HMM_FILE" "$temp_file" 2>> "$LOG_FILE"
    if [ $? -ne 0 ]; then
        echo "❌ Error: hmmsearch failed for $FASTA_FILE × $HMM_FILE" | tee -a "$LOG_FILE"
        rm -f "$tblout_file" "$domtblout_file"
    else
        echo "✅ Finished: $FASTA_FILE"
    fi

    rm -f "$temp_file"
}
export -f process_hmm_fasta

# 遍历每个 HMM，针对所有 FASTA 并行运行
find "$HMM_DIR" -name "*.hmm" | while read -r HMM_FILE; do
    find "$FASTA_DIR" -type f \( -name "*.faa" -o -name "*.faa.gz" -o -name "*.fasta" -o -name "*.fasta.gz" \) \
        | parallel -j "$PARALLEL_JOBS" process_hmm_fasta "$HMM_FILE" {} "$OUTPUT_DIR" "$HMM_THREADS"
done

echo "✅ All hmmsearch tasks completed. Results saved in: $OUTPUT_DIR"

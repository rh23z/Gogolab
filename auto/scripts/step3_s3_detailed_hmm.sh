#!/bin/bash

# 检查输入参数
if [ "$#" -lt 5 ]; then
    echo "Usage: $0 <hmm_directory> <fasta_directory> <output_directory> <threads> [Score_threshold]"
    echo "Example: $0 hmm_dir fasta_dir out_dir 20 6"
    exit 1
fi

# 获取输入参数
HMM_DIR="$1"
FASTA_DIR="$2"
OUTPUT_DIR="$3"
THREADS="$4"
SCORE_THRESHOLD="${5:-0}"  # 默认 score 阈值为 0

# 检查目录
[ ! -d "$HMM_DIR" ] && echo "Error: HMM directory '$HMM_DIR' not found." && exit 1
[ ! -d "$FASTA_DIR" ] && echo "Error: FASTA directory '$FASTA_DIR' not found." && exit 1

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 单文件处理函数
process_hmm_fasta() {
    HMM_FILE="$1"
    FASTA_FILE="$2"
    OUTPUT_DIR="$3"

    HMMNAME=$(basename "$HMM_FILE" .hmm)
    BASENAME=$(basename "$FASTA_FILE")
    BASENAME=${BASENAME%.fasta}

    OUTPUT_PATH="$OUTPUT_DIR/$HMMNAME"
    mkdir -p "$OUTPUT_PATH"
    LOG_FILE="$OUTPUT_PATH/hmmsearch_errors.log"

    tblout_file="$OUTPUT_PATH/${BASENAME}.tblout"
    domtblout_file="$OUTPUT_PATH/${BASENAME}.domtblout"
    temp_file="$OUTPUT_PATH/${BASENAME}_cleaned.faa"

    # 跳过已存在的结果
    if [ -s "$tblout_file" ] && [ -s "$domtblout_file" ]; then
        echo "Skipping $domtblout_file as results already exist."
        return
    fi

    # 清理 FASTA（去除非法字符 '-'）
    if [[ "$FASTA_FILE" == *.gz ]]; then
        zcat "$FASTA_FILE" | sed '/^>/! s/-//g' > "$temp_file"
    else
        sed '/^>/! s/-//g' "$FASTA_FILE" > "$temp_file"
    fi

    echo "Running hmmsearch: HMM=$HMM_FILE FASTA=$FASTA_FILE T=$SCORE_THRESHOLD"

    hmmsearch --cpu "$THREADS" \
        --domtblout "$domtblout_file" \
        --tblout "$tblout_file" \
        -T "$SCORE_THRESHOLD" \
        "$HMM_FILE" "$temp_file" 2>> "$LOG_FILE"

    if [ $? -ne 0 ]; then
        echo "Error: hmmsearch failed for '$FASTA_FILE' with model '$HMM_FILE'." | tee -a "$LOG_FILE"
        rm -f "$tblout_file" "$domtblout_file"
    else
        echo "✔ hmmsearch completed for $FASTA_FILE"
    fi

    rm -f "$temp_file"
}

export -f process_hmm_fasta
export SCORE_THRESHOLD THREADS

# 并行运行所有 hmm × fasta
find "$HMM_DIR" -name "*.hmm" | while read -r HMM_FILE; do
    find "$FASTA_DIR" -name "*.faa*" -o -name "*.fasta" | \
        parallel -j "$THREADS" process_hmm_fasta "$HMM_FILE" {} "$OUTPUT_DIR"
done

echo "✅ hmmsearch completed. Results in '$OUTPUT_DIR'."

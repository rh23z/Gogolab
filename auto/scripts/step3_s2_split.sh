#!/bin/bash
set -euo pipefail

# ========== 检查输入参数 ==========
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <input_dir> <output_dir> <split_parts> [merged_filename]"
    echo "Example: $0 /path/to/fasta /path/to/output 20 step3_merged.fasta"
    exit 1
fi

input_dir="$1"
output_dir="$2"
split_parts="$3"
merged_fasta="${4:-step3_merged.fasta}"  # 可选参数，默认为 step3_merged.fasta

# 检查目录
if [ ! -d "$input_dir" ]; then
    echo "[错误] 输入目录不存在: $input_dir"
    exit 1
fi

# 创建输出目录
mkdir -p "$output_dir"

# 设置合并文件路径 → 输入目录的母目录
parent_dir="$(dirname "$input_dir")"
merged_fasta_path="$parent_dir/$merged_fasta"

# 如果文件已存在，先清空
> "$merged_fasta_path"

# 合并所有 fasta 文件
echo "🔄 合并所有 .fasta 文件到: $merged_fasta_path"
find "$input_dir" -type f -name "*.fasta" -exec cat {} + >> "$merged_fasta_path"

# 检查合并结果
if [ ! -s "$merged_fasta_path" ]; then
    echo "[错误] 合并失败，未生成有效的 $merged_fasta_path"
    exit 1
fi

# 拆分
echo "✂️  拆分为 $split_parts 份，输出到: $output_dir"
seqkit split -p "$split_parts" -O "$output_dir" "$merged_fasta_path"

echo "✅ 完成"

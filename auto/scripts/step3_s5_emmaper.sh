#!/usr/bin/env bash
set -Eeuo pipefail

# 默认值
THREADS=32
JOBS=20
INPUT_DIR=""
OUTPUT_DIR=""
EMAPPER=""

# 打印帮助信息
print_usage() {
  cat <<USAGE
Usage: $(basename "$0") -i INPUT_DIR -o OUTPUT_DIR --emapper PATH [options]

Required:
  -i, --input DIR        输入包含 FASTA 的目录
  -o, --output DIR       输出目录（每个输入 FASTA 生成一组 emapper 结果）
  --emapper PATH         emapper.py 的路径

Optional:
  -t, --threads N        每个 emapper 进程使用的线程数（默认: ${THREADS})
  -j, --jobs N           同时运行的 emapper 任务数（默认: ${JOBS})
USAGE
}

# 参数解析
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)   INPUT_DIR="$2"; shift 2;;
    -o|--output)  OUTPUT_DIR="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    -j|--jobs)    JOBS="$2"; shift 2;;
    --emapper)    EMAPPER="$2"; shift 2;;
    -h|--help)    print_usage; exit 0;;
    *) echo "[ERR] Unknown arg: $1"; print_usage; exit 1;;
  esac
done

# 参数检查
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$EMAPPER" ]]; then
  echo "[ERR] 缺少必需参数 -i/-o/--emapper"
  print_usage
  exit 1
fi
[[ ! -d "$INPUT_DIR" ]] && echo "[ERR] 输入目录不存在: $INPUT_DIR" && exit 1
[[ ! -f "$EMAPPER" ]] && echo "[ERR] emapper.py 不存在: $EMAPPER" && exit 1

command -v python >/dev/null || { echo "[ERR] 找不到 python 解释器"; exit 1; }
command -v parallel >/dev/null || { echo "[ERR] 缺少 GNU parallel"; exit 1; }

mkdir -p "$OUTPUT_DIR"

# 收集所有 .fasta / .fa 文件
mapfile -t FILES < <(find "$INPUT_DIR" -type f \( -name "*.fasta" -o -name "*.fa" \) | sort)
[[ ${#FILES[@]} -eq 0 ]] && echo "[WARN] 未找到任何 fasta 文件" && exit 0

echo "[INFO] 共 ${#FILES[@]} 个文件 | 并发任务: $JOBS | 每任务线程: $THREADS"

# 运行函数
run_emapper() {
  local fasta_file="$1"
  local outdir="$2"
  local emapper_py="$3"
  local threads="$4"

  local base_name
  base_name=$(basename "$fasta_file")
  base_name=${base_name%.fasta}
  base_name=${base_name%.fa}

  local out_prefix="${outdir}/${base_name}"
  local log_file="${out_prefix}.emapper.log"

  echo "[RUN ] emapper -> $base_name"
  if ! $emapper_py -i "$fasta_file" -o "$out_prefix" \
       --itype proteins\
       -m diamond \
       --cpu "$threads" >"$log_file" 2>&1; then
    echo "[FAIL] $base_name，日志: $log_file"
    return 1
  fi
  echo "[DONE] $base_name"
}
export -f run_emapper
export OUTPUT_DIR THREADS EMAPPER

# 并行运行
parallel --halt now,fail=1 -j "$JOBS" \
  run_emapper {} "$OUTPUT_DIR" "$EMAPPER" "$THREADS" ::: "${FILES[@]}"

echo "✅ 所有 emapper 任务已完成"

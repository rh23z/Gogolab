#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integrate HMMER *domtblout* files (Step3 detailed HMM → integrate) — CLI version

功能：
- 扫描给定目录（可递归）下的所有 *.domtblout[.gz]* 文件；
- 并行解析，整合为一个 TSV；
- 列名采用 HMMER domtblout 标准 22 列 + description；
- 额外附带列：target_file（文件名或相对路径）。

用法示例：

python step3_s2_detailed_hmm_inte_cli.py \
  -i /work/zhangrh/couple_procject/result/cas9/flanking/Pfam-A \
  -o /work/zhangrh/couple_procject/result/cas9/flanking/step3_merged_pfam.tsv \
  --recursive --workers 128

可选：
  --ext .domtblout        # 默认
  --include-gz            # 同时匹配 .domtblout.gz
  --relative-path         # target_file 使用相对路径（相对于输入根目录）

"""
from __future__ import annotations
import os
import gzip
import argparse
from pathlib import Path
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Iterable, List

# HMMER domtblout 标准字段列表（22 + description）
COLNAMES = [
    "target_name", "target_accession", "tlen",
    "query_name", "query_accession", "qlen",
    "full_seq_Evalue", "full_seq_score", "full_seq_bias",
    "domain_num", "domain_of", "c_Evalue", "i_Evalue",
    "domain_score", "domain_bias",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
    "acc", "description"
]


def _open_text(path: Path):
    """Open text file, supporting optional .gz."""
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def parse_domtblout(domtblout_path: Path) -> pd.DataFrame | None:
    try:
        with _open_text(domtblout_path) as f:
            lines = [line.strip() for line in f if not line.startswith('#') and line.strip()]
        if not lines:
            return None

        records: List[List[str]] = []
        for line in lines:
            parts = line.split()
            if len(parts) < 22:
                continue
            fixed = parts[:22]
            description = ' '.join(parts[22:])
            fixed.append(description)
            records.append(fixed)

        if not records:
            return None

        df = pd.DataFrame(records, columns=COLNAMES)
        return df
    except Exception as e:
        print(f"❌ Error reading {domtblout_path}: {e}")
        return None


def _iter_domtblout_files(root: Path, recursive: bool, ext: str, include_gz: bool) -> Iterable[Path]:
    patterns = [f"*{ext}"]
    if include_gz:
        patterns.append(f"*{ext}.gz")
    if recursive:
        for pat in patterns:
            yield from root.rglob(pat)
    else:
        for pat in patterns:
            yield from root.glob(pat)


def main():
    ap = argparse.ArgumentParser(description='Integrate HMMER domtblout files into a single TSV (parallel).')
    ap.add_argument('-i', '--input-dir', required=True, help='输入目录，包含 *.domtblout[.gz] 文件')
    ap.add_argument('-o', '--output', required=True, help='输出 TSV 文件路径 (例如 step3_merged_pfam.tsv)')
    ap.add_argument('--ext', default='.domtblout', help='文件扩展名（默认 .domtblout）')
    ap.add_argument('--recursive', action='store_true', help='递归查找子目录')
    ap.add_argument('--include-gz', action='store_true', help='同时解析 .domtblout.gz')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 8, help='并行进程数（默认=CPU核数）')
    ap.add_argument('--relative-path', action='store_true', help='target_file 使用相对路径（相对于输入目录）')
    args = ap.parse_args()

    in_root = Path(args.input_dir).resolve()
    out_path = Path(args.output).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    domtblout_files = list(_iter_domtblout_files(in_root, args.recursive, args.ext, args.include_gz))
    if not domtblout_files:
        print(f"⚠️ No *{args.ext}[.gz] files found under: {in_root}")
        return

    print(f"🔎 Found {len(domtblout_files)} files. Parsing with workers={args.workers} …")
    merged_df_list: List[pd.DataFrame] = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(parse_domtblout, p): p for p in domtblout_files}
        for fut in as_completed(futures):
            p = futures[fut]
            df = fut.result()
            if df is None:
                continue
            # 附加 target_file 列
            df = df.copy()
            if args.relative_path:
                df['target_file'] = str(p.relative_to(in_root))
            else:
                df['target_file'] = p.name
            merged_df_list.append(df)

    if not merged_df_list:
        print("⚠️ No valid records to write.")
        return

    merged = pd.concat(merged_df_list, ignore_index=True)
    # 输出单个 TSV（避免原脚本中的双 .tsv）
    if not out_path.suffix:
        out_path = out_path.with_suffix('.tsv')
    merged.to_csv(out_path, sep='\t', index=False)
    print(f"✅ Written to: {out_path}")


if __name__ == '__main__':
    main()

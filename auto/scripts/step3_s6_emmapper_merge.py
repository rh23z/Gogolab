#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step3 S6 — Merge eggNOG-mapper results (CLI)

把一个目录下的 *.emapper.annotations 与 *.seed_orthologs 合并为单一 TSV。
- 支持递归、可选并行读取、可选 .gz（按扩展名自动判断）。
- annotations 与 seed_orthologs 均按 query_name 合并（left join）。

示例：
  python step3_s6_emapper_merge_cli.py \
    -i /work/zhangrh/couple_procject/result/cas9/flanking/egg \
    -o /work/zhangrh/couple_procject/result/cas9/flanking/step3_merged_egg.tsv \
    --recursive --workers 64
"""
from __future__ import annotations
import os
import argparse
from pathlib import Path
from typing import List, Optional
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

# ========== 期望列定义（与现有脚本保持一致） ==========
ANN_COLUMNS = [
    'query_name','seed_ortholog','evalue','score','eggNOG_OGs','max_annot_lvl','COG_category',
    'Description','Preferred_name','GOs','EC','KEGG_ko','KEGG_Pathway','KEGG_Module','KEGG_Reaction',
    'KEGG_rclass','BRITE','KEGG_TC','CAZy','BiGG_Reaction','PFAMs'
]
ORTH_COLUMNS = ['query_name','sseqid','evalue_orth','bitscore','qstart','qend','sstart','send','pident','qcov','scov']

# ========== 读文件工具 ==========

def _read_table_nohdr(path: Path) -> Optional[pd.DataFrame]:
    """读取去注释（#）的制表文件，无表头。空表返回 None。支持 .gz 由 pandas 自动处理。"""
    try:
        df = pd.read_csv(path, sep='\t', comment='#', header=None, dtype=str)
        if df.empty:
            return None
        return df
    except Exception as e:
        print(f"[WARN] 读取失败: {path} -> {e}")
        return None


def _gather_files(root: Path, pattern: str, recursive: bool) -> List[Path]:
    if recursive:
        return sorted(root.rglob(pattern))
    return sorted(root.glob(pattern))

# ========== 并行读取集合 ==========

def _batch_read(paths: List[Path], expected_cols: int, workers: int) -> List[pd.DataFrame]:
    if not paths:
        return []
    out: List[pd.DataFrame] = []
    if workers <= 1:
        for p in paths:
            df = _read_table_nohdr(p)
            if df is None:
                continue
            if df.shape[1] < expected_cols:
                print(f"[WARN] 列数({df.shape[1]})少于期望({expected_cols}): {p} -> 跳过")
                continue
            out.append(df.iloc[:, :expected_cols].copy())
        return out

    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_read_table_nohdr, p): p for p in paths}
        for fut in as_completed(futs):
            p = futs[fut]
            df = fut.result()
            if df is None:
                continue
            if df.shape[1] < expected_cols:
                print(f"[WARN] 列数({df.shape[1]})少于期望({expected_cols}): {p} -> 跳过")
                continue
            out.append(df.iloc[:, :expected_cols].copy())
    return out

# ========== 主流程 ==========

def main():
    ap = argparse.ArgumentParser(description='Merge eggNOG-mapper annotations and seed_orthologs to one TSV.')
    ap.add_argument('-i','--input', required=True, help='输入目录，包含 *.emapper.annotations / *.seed_orthologs')
    ap.add_argument('-o','--output', required=True, help='输出 TSV 路径')
    ap.add_argument('--recursive', action='store_true', help='递归搜索子目录')
    ap.add_argument('--workers', type=int, default=1, help='并行读取进程数（默认 1，建议 8~64）')
    args = ap.parse_args()

    in_dir = Path(args.input).resolve()
    out_path = Path(args.output).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    ann_paths = _gather_files(in_dir, '*.emapper.annotations*', args.recursive)
    orth_paths = _gather_files(in_dir, '*.seed_orthologs*', args.recursive)

    if not ann_paths:
        print(f"[ERR] 找不到 annotations：{in_dir}")
        return
    if not orth_paths:
        print(f"[WARN] 找不到 seed_orthologs：{in_dir}（将只输出 annotations）")

    print(f"[INFO] annotations: {len(ann_paths)} | seed_orthologs: {len(orth_paths)} | workers={args.workers}")

    ann_dfs = _batch_read(ann_paths, expected_cols=len(ANN_COLUMNS), workers=args.workers)
    if not ann_dfs:
        print('[ERR] 无有效 annotations 记录，退出。')
        return
    ann_merged = pd.concat(ann_dfs, ignore_index=True)
    ann_merged.columns = ANN_COLUMNS

    if orth_paths:
        orth_dfs = _batch_read(orth_paths, expected_cols=len(ORTH_COLUMNS), workers=args.workers)
        if orth_dfs:
            orth_merged = pd.concat(orth_dfs, ignore_index=True)
            orth_merged.columns = ORTH_COLUMNS
            final_df = pd.merge(ann_merged, orth_merged, on='query_name', how='left')
        else:
            print('[WARN] 无有效 seed_orthologs 记录，跳过合并。')
            final_df = ann_merged
    else:
        final_df = ann_merged

    final_df.to_csv(out_path, sep='\t', index=False)
    print(f"[OK] 写出：{out_path}")


if __name__ == '__main__':
    main()

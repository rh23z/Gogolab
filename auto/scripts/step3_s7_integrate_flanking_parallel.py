#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step3 S7 — Integrate flanking annotations (CLI, multi-processing accelerated)

在 s7 CLI 的基础上，加入 --workers 并行展开 flanking_segments（分块并行），
随后与 PFAM / eggNOG 合并并写出。

示例：
  python step3_s7_integrate_flanking_cli_mt.py \
    --pfam /work/zhangrh/couple_procject/result/cas9/flanking/step3_merged_pfam.tsv \
    --emapper /work/zhangrh/couple_procject/result/cas9/flanking/step3_merged_egg.tsv \
    --step3-df /work/zhangrh/couple_procject/result/cas9/processed_df/step3_df.tsv \
    -o /work/zhangrh/couple_procject/result/cas9/processed_df/step3_df_position.tsv \
    --workers 64
"""
from __future__ import annotations
import argparse
import ast
from pathlib import Path
from typing import Iterable, List, Dict, Any
from concurrent.futures import ProcessPoolExecutor, as_completed
import math

import pandas as pd


# ------------------------------
# I/O helpers
# ------------------------------

def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Not found: {path}")
    # pandas 会根据后缀自动识别压缩；engine='c' 更快
    return pd.read_table(path, engine='c')


# ------------------------------
# Load PFAM (domtblout merged)
# ------------------------------

def load_pfam_grouped(pfam_path: Path) -> pd.DataFrame:
    p_df = pd.read_table(pfam_path)

    # 重命名列
    p_df = p_df.rename(columns={
        'query_name': 'pfam_query_name',
        'target_name': 'flanking_id',
        'full_seq_Evalue': 'pfam_evalue',
        'full_seq_score': 'pfam_score'
    })

    # 保留需要的列并按 flanking_id 分组，将其余列变为 list
    grouped_df = (
        p_df[['flanking_id', 'pfam_evalue', 'pfam_score', 'pfam_query_name']]
        .groupby('flanking_id', as_index=False)
        .agg(list)
    )

    return grouped_df



# ------------------------------
# Load eggNOG-mapper merged
# ------------------------------

def load_emapper(emapper_path: Path) -> pd.DataFrame:
    em_df = _read_tsv(emapper_path)
    expect = ['query_name', 'evalue', 'score', 'Description', 'PFAMs',
              'sseqid', 'qstart', 'qend', 'sstart', 'send', 'pident', 'qcov', 'scov']
    for c in expect:
        if c not in em_df.columns:
            em_df[c] = pd.NA
    em_df = em_df[expect]
    em_df = em_df.rename(columns={
        'query_name': 'flanking_id',
        'evalue': 'emapper_evalue',
        'qstart': 'emapper_protein_start',
        'qend': 'emapper_protein_end',
        'qcov': 'emapper_protein_cov',
        'PFAMs':'emmapper_PFAMs',
        'Description': 'emmapper_Description',
    })
    keep = ['flanking_id', 'emapper_evalue', 'score', 'emmapper_Description', 'emmapper_PFAMs',
            'sseqid', 'emapper_protein_start', 'emapper_protein_end',
            'sstart', 'send', 'pident', 'emapper_protein_cov', 'scov']
    # 去重以防 emapper 多次输出导致同一个 flanking_id 多行
    return em_df[keep].drop_duplicates(subset=['flanking_id']).copy()


# ------------------------------
# Expand step3_df flanking segments — parallel
# ------------------------------

def _safe_parse_segments(x: Any) -> List:
    if isinstance(x, list):
        return x
    if pd.isna(x):
        return []
    if isinstance(x, str):
        s = x.strip()
        if not s:
            return []
        try:
            return ast.literal_eval(s)
        except Exception:
            return []
    return []


def _expand_rows(rows: Iterable[Dict[str, Any]], keep_cols: List[str]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for row in rows:
        segs = row.get('flanking_segments', []) or []
        for item in segs:
            if not isinstance(item, (list, tuple)) or len(item) < 4:
                continue
            out.append({
                **{c: row.get(c) for c in keep_cols},
                'flanking_id': item[0],
                'flanking_start': int(item[1]),
                'flanking_end': int(item[2]),
                'flanking_strand': int(item[3]),
            })
    return out


def expand_step3_parallel(step3_df: pd.DataFrame, keep_cols: List[str], workers: int = 1, chunksize: int = 5000) -> pd.DataFrame:
    step3_df = step3_df.copy()
    step3_df['flanking_segments'] = step3_df['flanking_segments'].apply(_safe_parse_segments)

    if workers <= 1 or len(step3_df) <= chunksize:
        rows = step3_df.to_dict('records')
        expanded = _expand_rows(rows, keep_cols)
        return pd.DataFrame(expanded)

    # 分块
    n = len(step3_df)
    n_chunks = math.ceil(n / chunksize)
    futures = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for i in range(n_chunks):
            start = i * chunksize
            end = min((i + 1) * chunksize, n)
            chunk_rows = step3_df.iloc[start:end].to_dict('records')
            futures.append(ex.submit(_expand_rows, chunk_rows, keep_cols))
        expanded_parts: List[List[Dict[str, Any]]] = []
        for fut in as_completed(futures):
            expanded_parts.append(fut.result())
    # 合并
    flat: List[Dict[str, Any]] = [r for part in expanded_parts for r in part]
    return pd.DataFrame(flat)


# ------------------------------
# Distance calculation
# ------------------------------

def add_relative_distance(df_expanded: pd.DataFrame) -> pd.DataFrame:
    # 向量化计算，避免逐行 apply
    ps = pd.to_numeric(df_expanded['prot_start'], errors='coerce')
    pe = pd.to_numeric(df_expanded['prot_end'], errors='coerce')
    fs = pd.to_numeric(df_expanded['flanking_start'], errors='coerce')
    fe = pd.to_numeric(df_expanded['flanking_end'], errors='coerce')

    # 默认 0（重叠）
    rel = pd.Series(0, index=df_expanded.index, dtype='int64')
    # 上游：fe < ps -> fe - ps
    mask_up = (fe < ps)
    rel.loc[mask_up] = (fe - ps)[mask_up]
    # 下游：fs > pe -> fs - pe
    mask_down = (fs > pe)
    rel.loc[mask_down] = (fs - pe)[mask_down]

    df_expanded['relative_distance'] = rel
    return df_expanded


# ------------------------------
# CLI
# ------------------------------

def main():
    ap = argparse.ArgumentParser(description='Integrate PFAM & eggNOG-mapper annotations onto flanking segments table (multi-processing).')
    ap.add_argument('--pfam', required=True, help='PFAM merged TSV (from domtblout integration)')
    ap.add_argument('--emapper', required=True, help='eggNOG-mapper merged TSV (from S6)')
    ap.add_argument('--step3-df', required=True, help='Step3 dataframe TSV with flanking_segments column')
    ap.add_argument('-o', '--out', required=True, help='Output TSV path')
    ap.add_argument('--keep-cols', default='target_name,target_file,prot_start,prot_end,source',
                    help='Comma-separated columns to keep from step3_df (default matches original)')
    ap.add_argument('--workers', type=int, default=32, help='并行进程数（>1 开启分块并行）')
    ap.add_argument('--chunksize', type=int, default=5000, help='每个进程处理的行数块大小（默认 5000）')
    args = ap.parse_args()

    pfam_path = Path(args.pfam)
    emapper_path = Path(args.emapper)
    step3_path = Path(args.step3_df)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    keep_cols = [c.strip() for c in args.keep_cols.split(',') if c.strip()]

    # 1) 读 & 重命名
    p_df = load_pfam_grouped(pfam_path)
    e_df = load_emapper(emapper_path)

    # 2) 读 step3 & 并行展开
    step3_df = _read_tsv(step3_path)
    if 'flanking_segments' not in step3_df.columns:
        raise ValueError("Input step3_df missing column: flanking_segments")
    expanded = expand_step3_parallel(step3_df, keep_cols, workers=args.workers, chunksize=args.chunksize)

    if expanded.empty:
        print('[WARN] 展开结果为空，仍将写出含表头的空文件。')
        expanded.to_csv(out_path, sep='\t', index=False)
        return

    # 3) 相对距离
    expanded = add_relative_distance(expanded)

    # 4) 合并 emapper / pfam（按 flanking_id 左连接）
    out = expanded.merge(
        e_df[['flanking_id', 'emapper_evalue', 'emapper_protein_start', 'emapper_protein_end', 'emapper_protein_cov','emmapper_PFAMs', 'emmapper_Description']],
        on='flanking_id', how='left'
    ).merge(
        p_df[['flanking_id', 'pfam_evalue', 'pfam_score', 'pfam_query_name']],
        on='flanking_id', how='left'
    )

    # 5) 写出
    if not out_path.suffix:
        out_path = out_path.with_suffix('.tsv')
    out.to_csv(out_path, sep='\t', index=False)
    print(f"[OK] Written: {out_path}")


if __name__ == '__main__':
    main()

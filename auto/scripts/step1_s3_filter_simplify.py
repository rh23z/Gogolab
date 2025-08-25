#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import argparse
import pandas as pd


def get_position(description):
    match = re.findall(r"#\s*(\d+)\s*#\s*(\d+)\s*#\s*(-?\d+)\s*#", description)
    if match:
        return tuple(map(int, match[0]))
    else:
        return (None, None, None)


def process_single_file(path, score_cutoff):
    df = pd.read_table(path)
    df = df[df['full_seq_score'].astype(float) > score_cutoff].copy()
    df['source'] = os.path.basename(path).replace('_domtblout_merged.tsv', '')
    df[['prot_start', 'prot_end', 'strand']] = df['description'].apply(lambda x: pd.Series(get_position(x)))

    filtered_out = path.replace('_domtblout_merged.tsv', '_domtblout_merged_filtered.tsv')
    df.to_csv(filtered_out, sep='\t', index=False)

    results = []
    for id_, group in df.groupby("target_name"):
        result = {
            'target_name': id_,
            'query_names': group['query_name'].tolist(),
            'full_seq_Evalue': group['full_seq_Evalue'].tolist(),
            'full_seq_score': group['full_seq_score'].tolist(),
            'full_seq_bias': group['full_seq_bias'].tolist(),
            'i_Evalue': group['i_Evalue'].tolist(),
            'domain_score': group['domain_score'].tolist(),
            'domain_bias': group['domain_bias'].tolist(),
            'domain_mdl_len': group['qlen'].tolist(),
            'target_len': group['tlen'].iloc[0],
            'target_file': group['target_file'].iloc[0] if 'target_file' in group else None,
            'prot_start': group['prot_start'].iloc[0],
            'prot_end': group['prot_end'].iloc[0],
            'strand': group['strand'].iloc[0],
            'domain_alignment_positions': list(zip(group['env_from'], group['env_to'])),
            'domain_mdl_hit_coverage': ((group['env_to'] - group['env_from']) / group['qlen']).tolist(),
            'source': group['source'].iloc[0]
        }
        results.append(result)

    summary_out = path.replace('_domtblout_merged.tsv', '_summary.tsv')
    pd.DataFrame(results).to_csv(summary_out, sep='\t', index=False)
    print(f"✅ Processed {os.path.basename(path)} → {os.path.basename(summary_out)}")


def main(args):
    input_dir = args.input_dir
    output_prefix = args.output_prefix
    score_cutoff = args.score_cutoff

    # 处理每个 *_merged.tsv 文件
    for file in os.listdir(input_dir):
        if file.endswith('_domtblout_merged.tsv'):
            full_path = os.path.join(input_dir, file)
            process_single_file(full_path, score_cutoff)

    # 合并所有 *_summary.tsv
    combined = pd.DataFrame()
    for file in os.listdir(input_dir):
        if file.endswith('_summary.tsv'):
            df = pd.read_table(os.path.join(input_dir, file))
            combined = pd.concat([combined, df], ignore_index=True)

    final_out = f"{output_prefix}.tsv"
    combined.to_csv(final_out, sep='\t', index=False)
    print(f"✅ Combined summary saved to: {final_out}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Step1.3: Filter and simplify domtblout results')
    parser.add_argument('--input_dir', required=True, help='Directory containing *_domtblout_merged.tsv files')
    parser.add_argument('--output_prefix', default='step1_filtered_df', help='Prefix for final combined TSV output (default: step1_filtered_df)')
    parser.add_argument('--score_cutoff', type=float, default=16.0, help='Minimum full_seq_score (default: 16.0)')
    args = parser.parse_args()
    main(args)

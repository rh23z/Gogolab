#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed


# ========== 解析 domtblout 文件 ==========
def parse_domtblout(domtblout_path):
    try:
        with open(domtblout_path, 'r') as f:
            lines = [line.strip() for line in f if not line.startswith('#') and line.strip()]
        if not lines:
            return None

        records = []
        for line in lines:
            parts = line.split()
            if len(parts) < 22:
                continue
            fixed = parts[:22]
            description = ' '.join(parts[22:])
            fixed.append(description)
            records.append(fixed)

        colnames = [
            "target_name", "target_accession", "tlen",
            "query_name", "query_accession", "qlen",
            "full_seq_Evalue", "full_seq_score", "full_seq_bias",
            "domain_num", "domain_of", "c_Evalue", "i_Evalue",
            "domain_score", "domain_bias",
            "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
            "acc", "description"
        ]
        df = pd.DataFrame(records, columns=colnames)
        df['target_file'] = os.path.basename(domtblout_path)
        return df

    except Exception as e:
        print(f"❌ Error reading {domtblout_path}: {e}")
        return None


# ========== 主函数 ==========
def main(args):
    input_dir = args.input_dir
    output_dir = args.output_dir
    threads = args.threads

    os.makedirs(output_dir, exist_ok=True)

    for results_dir in os.listdir(input_dir):
        full_results_dir = os.path.join(input_dir, results_dir)
        if not os.path.isdir(full_results_dir):
            continue

        print(f"🔍 Scanning directory: {full_results_dir}")
        domtblout_files = []
        for root, _, files in os.walk(full_results_dir):
            for file in files:
                if file.endswith('.domtblout'):
                    domtblout_files.append(os.path.join(root, file))

        print(f"📄 Found {len(domtblout_files)} .domtblout files")

        if not domtblout_files:
            print(f"⚠️ No valid .domtblout files in {full_results_dir}")
            continue

        merged_df_list = []

        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {executor.submit(parse_domtblout, path): path for path in domtblout_files}
            for future in as_completed(futures):
                df = future.result()
                if df is not None:
                    merged_df_list.append(df)

        if merged_df_list:
            merged_df = pd.concat(merged_df_list, ignore_index=True)
            output_file = os.path.join(output_dir, f"{results_dir}_domtblout_merged.tsv")
            merged_df.to_csv(output_file, sep='\t', index=False)
            print(f"✅ Written to: {output_file}")
        else:
            print(f"⚠️ No valid records to write for {results_dir}")


# ========== 命令行入口 ==========
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step1.2: Merge multiple .domtblout files in parallel")
    parser.add_argument('--input_dir', required=True, help='Parent directory containing subfolders with domtblout files')
    parser.add_argument('--output_dir', required=True, help='Directory to save merged TSV files')
    parser.add_argument('--threads', type=int, default=64, help='Number of parallel threads (default: 64)')
    args = parser.parse_args()

    main(args)

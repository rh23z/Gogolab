#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
命令行工具：根据 summary.tsv 中的 source/file_name/seq_id 信息，
从不同目录的 fasta|faa 文件里提取目标序列，最终合并输出一个 fasta。
"""

import argparse
import sys
import os
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd

# source 与 fasta 文件根目录的映射
SOURCE_DICT = {
    'archaea_assembly_summary': "/work/data_share/NCBI_zrh/archaea_assembly_summary_prodigal/",
    'bacteria_assembly_summary': "/work/data_share/NCBI_zrh/bacteria_assembly_summary_prodigal/",
    'BGItrans': "/work/data/BGItrans/final_bins_faa/",
    'metagenomes_assembly_summary': "/work/data_share/NCBI_zrh/metagenomes_assembly_summary_prodigal/",
    'phagescope': '/work/data/phage/phagescope_all_prodigal/',
    'viral_assembly_summary': "/work/data_share/NCBI_zrh/viral_assembly_summary_prodigal/",
    'IMGM': '/work/data/IMGM_1025/faa/',
    'MGnify': '/work/data/MGnify/faa/',
    'CRBC': '/work/data/CRBC_genome/faa/',
    'IMGM_metagenome': '/work/data/IMGM_metagenome/faa/'
}


def extract_sequences(source, file_name, need_ids, log_file):
    root = SOURCE_DICT.get(source)
    if root is None:
        print(f"[WARN] 未知 source: {source}", file=sys.stderr)
        return []

    fasta_path = os.path.join(root, file_name)
    if not os.path.exists(fasta_path):
        print(f"[WARN] 文件不存在: {fasta_path}", file=sys.stderr)
        return []

    out_records = []
    not_found = []

    for rec in SeqIO.parse(fasta_path, "fasta"):
        if rec.id in need_ids:
            out_records.append(rec)
            need_ids.remove(rec.id)

    not_found.extend(need_ids)
    if not_found:
        with open(log_file, 'a') as log:
            for seq_id in not_found:
                log.write(f"{source}\t{file_name}\t{seq_id}\n")

    return out_records


def main(args):
    input_tsv = args.input
    output_prefix = args.output_prefix
    log_dir = args.log_dir

    log_file = os.path.join(log_dir, "step2_s1_not_found.log")
    output_fasta = output_prefix+ ".fasta"

    print(f"🧾 读取输入表：{input_tsv}")
    df = pd.read_table(input_tsv)
    df['target_file'] = df['target_file'].str.replace(".domtblout", ".faa", regex=False)
    need_dict = df.groupby(["source", "target_file"])["target_name"].apply(set).to_dict()

    print(f"🚀 启动并行提取（线程数: {args.threads}）")
    all_records = []
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(extract_sequences, source, file_name, ids.copy(), log_file)
                   for (source, file_name), ids in need_dict.items()]
        for future in as_completed(futures):
            all_records.extend(future.result())

    print(f"🧬 提取完毕，共提取序列：{len(all_records):,} 条")
    SeqIO.write(all_records, output_fasta, "fasta")
    print(f"✅ 输出写入：{os.path.abspath(output_fasta)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="从多个来源的 FASTA 文件中提取目标序列并输出合并 FASTA")
    parser.add_argument('-i', '--input', required=True, help="包含 source、file、seq_id 的 tsv 文件")
    parser.add_argument('-o', '--output_prefix', required=True, help="输出目录路径（用于保存 .fasta 和日志）")
    parser.add_argument('-t', '--threads', type=int, default=32, help="并行线程数（默认: 32）")
    parser.add_argument('--log_dir', default='.', help="日志文件目录（默认: 当前目录）")
    args = parser.parse_args()
    main(args)

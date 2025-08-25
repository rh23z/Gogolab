#!/usr/bin/env python3
import argparse
import pandas as pd
import os
from Bio import SeqIO


def get_repseq_for_round(df, cluster_file, round_num, identity_thresholds):
    if round_num == 0:
        clustering_seqid = 'target_name'
    else:
        clustering_seqid = f'cluster_{identity_thresholds[round_num - 1]}_repseq'

    cluster_df = pd.read_table(cluster_file, header=None, names=['repseq_id', 'seqid'])
    df_merged = pd.merge(df, cluster_df, left_on=clustering_seqid, right_on='seqid', how='left')
    df_merged.rename(columns={'repseq_id': f'cluster_{identity_thresholds[round_num]}_repseq'}, inplace=True)
    df_merged.drop(columns=['seqid'], inplace=True)
    return df_merged


def process_clusters(input_dir, df, identity_thresholds):
    for round_num, identity in enumerate(identity_thresholds):
        cluster_file = os.path.join(input_dir, f"round{round_num + 1}_{identity}_cluster.tsv")
        if not os.path.isfile(cluster_file):
            raise FileNotFoundError(f"[错误] 找不到聚类结果文件: {cluster_file}")
        df = get_repseq_for_round(df, cluster_file, round_num, identity_thresholds)
    return df


def add_sequence_column(df, fasta_path):
    print(f"[INFO] 正在从 {fasta_path} 加载序列...")
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    missing = 0
    seqs = []
    for name in df["target_name"]:
        if name in seq_dict:
            seqs.append(str(seq_dict[name].seq))
        else:
            seqs.append("")
            missing += 1
    if missing > 0:
        print(f"[WARN] 有 {missing} 个 target_name 在 fasta 中未找到序列")
    df["seq"] = seqs
    return df


def main():
    parser = argparse.ArgumentParser(description="根据MMseqs聚类结果为每条记录添加代表序列列，并可提取FASTA序列")
    parser.add_argument("--input_df", required=True, help="输入的原始DF表 (包含 target_name 列)")
    parser.add_argument("--cluster_dir", required=True, help="聚类结果目录")
    parser.add_argument("--identity_thresholds", required=True, nargs='+', type=int,
                        help="identity 阈值列表，例如：100 98 90 70 50 30")
    parser.add_argument("--select_identity", required=True, nargs='+', type=int,
                        help="要提取代表序列的 identity 阈值列表，例如：50 70")
    parser.add_argument("--output_prefix", required=True, help="输出文件前缀")
    parser.add_argument("--input_fasta", required=False, help="可选：FASTA 文件用于提取 target_name 对应的序列")
    args = parser.parse_args()

    df = pd.read_table(args.input_df)
    df = process_clusters(args.cluster_dir, df, args.identity_thresholds)

    # 如提供 fasta，提取对应序列列
    if args.input_fasta:
        df = add_sequence_column(df, args.input_fasta)

    # 输出全表
    full_output = f"{args.output_prefix}_result.tsv"
    df.to_csv(full_output, sep='\t', index=False)

    # 输出指定 identity 下的代表序列子集
    for identity in args.select_identity:
        col_name = f'cluster_{identity}_repseq'
        if col_name not in df.columns:
            raise ValueError(f"[错误] 未找到 identity = {identity} 的代表序列列")
        sub_df = df[df['target_name'] == df[col_name]]
        sub_output = f"{args.output_prefix}_identity{identity}_result.tsv"
        sub_df.to_csv(sub_output, sep='\t', index=False)


if __name__ == "__main__":
    main()

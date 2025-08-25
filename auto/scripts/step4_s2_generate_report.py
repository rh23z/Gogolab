#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import math

def split_dataframe(df, num_parts):
    """将 DataFrame 平均拆分为 num_parts 个子 DataFrame"""
    chunk_size = math.ceil(len(df) / num_parts)
    return [df.iloc[i * chunk_size: (i + 1) * chunk_size] for i in range(num_parts)]

def main():
    parser = argparse.ArgumentParser(description="将一个 DataFrame 平均拆分为多个部分")
    parser.add_argument("--input", required=True, help="输入的 TSV 文件路径")
    parser.add_argument("--output_dir", required=True, help="输出子文件的目录")
    parser.add_argument("--num_parts", type=int, required=True, help="拆分的子文件数目")
    parser.add_argument("--prefix", default="split_part", help="输出文件名前缀 (默认: split_part)")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = pd.read_table(args.input)

    chunks = split_dataframe(df, args.num_parts)

    for idx, chunk in enumerate(chunks):
        output_path = os.path.join(args.output_dir, f"{args.prefix}_{idx+1}.tsv")
        chunk.to_csv(output_path, sep='\t', index=False)
        print(f"✅ Saved: {output_path} ({len(chunk)} rows)")

if __name__ == "__main__":
    main()

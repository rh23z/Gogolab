#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import ast
import argparse
import pandas as pd
import ast
from concurrent.futures import ProcessPoolExecutor


def main(args):
    df = pd.read_table(args.input_file)

    df = df[df['tlen'] > args.min_target_len].copy()

    filter_df1 = df[df['query_names'].apply(lambda x: all(keyword in x for keyword in args.and_filters))]
    filter_df1_list = filter_df1['target_name'].tolist()

    filter_df2 = df[df['query_names'].apply(lambda x: any(keyword in x for keyword in args.any_filters))]
    filter_df2_list = filter_df2['target_name'].tolist()

    print(f"✅ Filtered {len(filter_df1_list)} targets with AND filters")
    print(f"✅ Filtered {len(filter_df2_list)} targets with ANY filters")

    filter_df = pd.concat([filter_df1, filter_df2]).drop_duplicates()
    filter_df['query_names'] = filter_df['query_names'].apply(ast.literal_eval)
    filter_df['domain_mdl_hit_coverage'] = filter_df['domain_mdl_hit_coverage'].apply(ast.literal_eval)

    coverage_check_ids = [id_ for id_ in filter_df2_list if id_ not in filter_df1_list]

    def coverage_filter(row, conditions, threshold):

        for condition in conditions:
            if condition in row['query_names']:
                indices = [i for i, x in enumerate(row['query_names']) if x == condition]
                coverage_values = [row['domain_mdl_hit_coverage'][i] for i in indices]
                if any(cov > threshold for cov in coverage_values):
                    return True
        return False

    def parallel_filter(df, coverage_check_ids, conditions, threshold, threads):
        keep_rows = df[~df['target_name'].isin(coverage_check_ids)].index.tolist()

        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {
                executor.submit(coverage_filter, row, conditions, threshold): idx
                for idx, row in df[df['target_name'].isin(coverage_check_ids)].iterrows()
            }

            for future in tqdm(as_completed(futures), total=len(futures), desc="Processing rows"):
                idx = futures[future]
                if future.result():
                    keep_rows.append(idx)

        return df.loc[keep_rows]

    filter_df = parallel_filter(filter_df, coverage_check_ids, args.any_filters, args.cov_threshold, args.threads)



    filter_df.to_csv(args.output_file, sep='\t', index=False)
    print(f"✅ Filtered result saved to {args.output_file}")


if __name__ == '__main__':
    def parse_args():
        parser = argparse.ArgumentParser(
            description="Step1.4: Functional filter with AND/ANY domain logic + coverage + length"
        )
        parser.add_argument('--input_file', required=True, help='Input TSV file')
        parser.add_argument('--output_file', required=True, help='Output TSV file')
        parser.add_argument('--cov_threshold', type=float, default=0.2)
        parser.add_argument('--min_target_len', type=int, default=400)

        # 允许传入多个 AND 和 ANY filter
        parser.add_argument('--and_filters', nargs='*', default=[], help='List of filters that must be present in query_names')
        parser.add_argument('--any_filters', nargs='*', default=[], help='List of filters that should be present in query_names')

        parser.add_argument('--threads', type=int, default=4)

        return parser.parse_args()

    args = parse_args()
    main(args)


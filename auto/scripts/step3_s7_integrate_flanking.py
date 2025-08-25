import pandas as pd
import ast

# 读取 pfam 结果并重命名
p_df = pd.read_table('/work/zhangrh/couple_procject/result/cas9/flanking/step3_merged_pfam.tsv')
p_df.rename(columns={'query_name': 'pfam_query_name',
                     'target_name':'flanking_id',
                     'full_seq_Evalue':'pfam_evalue',
                     'full_seq_score':'pfam_score'}, inplace=True)

# 读取 emapper 结果并重命名列（保留 PFAMs 列）
em_df = pd.read_table('/work/zhangrh/couple_procject/result/cas9/flanking/step3_merged_egg.tsv')
em_df = em_df[['query_name', 'evalue', 'score', 'Description', 'PFAMs',
               'sseqid', 'qstart', 'qend', 'sstart', 'send', 'pident', 'qcov', 'scov']]
em_df.rename(columns={
    'query_name': 'flanking_id',
    'evalue': 'emapper_evalue',
    'qstart': 'emapper_protein_start',
    'qend': 'emapper_protein_end',
    'qcov': 'emapper_protein_cov'
}, inplace=True)

# 读取主 df
df = pd.read_table('/work/zhangrh/couple_procject/result/cas9/processed_df/step3_df.tsv')
df['flanking_segments'] = df['flanking_segments'].apply(
    lambda x: ast.literal_eval(x) if isinstance(x, str) else x
)
# 需要保留的原始列
keep_cols = [
    'target_name', 'target_file', 'prot_start', 'prot_end',
    'cluster_100_repseq', 'source'
]
# 拆分 flanking_segments
rows = []
for index, row in df.iterrows():
    flanking_segments = row['flanking_segments']
    for flanking_segments_item in flanking_segments:
        flanking_id = flanking_segments_item[0]

        if isinstance(flanking_segments_item, (list, tuple)) and len(flanking_segments_item) >= 4:
            rows.append({
                # 保留指定列
                **{col: row[col] for col in keep_cols if col in row},
                # 新拆分列
                'flanking_id': flanking_segments_item[0],
                'flanking_start': flanking_segments_item[1],
                'flanking_end': flanking_segments_item[2],
                'flanking_strand': flanking_segments_item[3]
            })

df_expanded = pd.DataFrame(rows)

# 计算具体距离
def calc_distance(row):
    if row['flanking_end'] < row['prot_start']:  # 上游
        return row['flanking_end'] - row['prot_start']
    elif row['flanking_start'] > row['prot_end']:  # 下游
        return row['flanking_start'] - row['prot_end']
    else:  # 重叠
        return 0

df_expanded['relative_distance'] = df_expanded.apply(calc_distance, axis=1)
df_expanded = pd.merge(
    df_expanded, em_df[['flanking_id', 'emapper_evalue', 'emapper_protein_start',
                        'emapper_protein_end', 'emapper_protein_cov']],
    left_on='flanking_id', right_on='flanking_id', how='left'
)
df_expanded = pd.merge(
    df_expanded, p_df[['flanking_id', 'pfam_evalue', 'pfam_score', 'pfam_query_name']],
    left_on='flanking_id', right_on='flanking_id', how='left'
)

# 保存展开结果
df_expanded.to_csv(
    '/work/zhangrh/couple_procject/result/cas9/processed_df/step3_df_posistion.tsv',
    sep='\t', index=False
)

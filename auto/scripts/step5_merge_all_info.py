import re
import pandas as pd
import os

def get_repseq_for_round(df, cluster_file, round_num, identity_thresholds):
    # 读取当前轮次的_cluster.tsv文件
    if round_num == 0:
        clustering_seqid = 'target_name'
    else:
        clustering_seqid = f'cluster_{identity_thresholds[round_num-1]}_repseq'
    cluster_df = pd.read_table(cluster_file, header=None, names=['repseq_id', 'seqid'])

    df_merged = pd.merge(df, cluster_df, left_on=clustering_seqid, right_on='seqid', how='left')

    df_merged.rename(columns={'repseq_id': f'cluster_{identity_thresholds[round_num]}_repseq'}, inplace=True)
    df_merged.drop(columns=['seqid'], inplace=True)
    return df_merged

def process_clusters(input_dir, df, identity_thresholds):
    for round_num, identity in enumerate(identity_thresholds):
        # 使用 os.path.join 来构造文件路径
        cluster_file = os.path.join(input_dir, f"round{round_num+1}_{identity}_cluster.tsv")
        df = get_repseq_for_round(df, cluster_file, round_num, identity_thresholds)

    return df

def parse_crispr_file(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    organism = None
    crispr_list = []
    current_crispr = None
    inside_table = False
    table_lines = []
    repeat_info = None  # 保存统计行数据

    for line in lines:
        line = line.rstrip('\n')

        # ORGANISM
        m_org = re.match(r'^ORGANISM:\s+(\S+)', line)
        if m_org:
            organism = m_org.group(1)
            continue

        # CRISPR 开头
        m_crispr = re.match(r'^CRISPR\s+(\d+)\s+Range:\s*(\d+)\s*-\s*(\d+)', line)
        if m_crispr:
            if current_crispr and table_lines:
                crispr_list.extend(parse_crispr_table(current_crispr, table_lines, organism, repeat_info))
                table_lines = []
                repeat_info = None

            current_crispr = {
                'crispr_num': int(m_crispr.group(1)),
                'start': int(m_crispr.group(2)),
                'end': int(m_crispr.group(3)),
            }
            inside_table = False
            continue

        if re.match(r'^POSITION\s+REPEAT\s+SPACER', line):
            inside_table = True
            continue

        if inside_table:
            if re.match(r'^-+', line):
                continue
            # Repeats: 6	Average Length: 35		Average Length: 28
            if line.startswith('Repeats:'):
                m_stats = re.search(r'Repeats:\s*(\d+)\s+Average Length:\s*(\d+)\s+Average Length:\s*(\d+)', line.replace('\t', ' '))
                if m_stats:
                    repeat_info = {
                        'repeat_count': int(m_stats.group(1)),
                        'repeat_avg_len': int(m_stats.group(2)),
                        'spacer_avg_len': int(m_stats.group(3)),
                    }
                continue
            table_lines.append(line)

    # 文件末尾最后一个CRISPR区块处理
    if current_crispr and table_lines:
        crispr_list.extend(parse_crispr_table(current_crispr, table_lines, organism, repeat_info))

    return pd.DataFrame(crispr_list)

def parse_crispr_table(crispr_info, lines, organism, repeat_info=None):
    records = []
    for line in lines:
        parts = re.split(r'\s{2,}|\t', line.strip())
        if len(parts) < 2:
            continue
        position = parts[0]
        repeat = parts[1]
        spacer = parts[2] if len(parts) > 2 else ''
        spacer_lengths = None
        if len(parts) > 3:
            m = re.findall(r'\d+', parts[3])
            spacer_lengths = list(map(int, m)) if m else None

        record = {
            'genome': organism,
            'crispr_num': crispr_info['crispr_num'],
            'crispr_start': crispr_info['start'],
            'crispr_end': crispr_info['end'],
            'position': int(position),
            'repeat_seq': repeat,
            'spacer_seq': spacer,
            'spacer_lengths': spacer_lengths,
            'repeat_count': repeat_info['repeat_count'] if repeat_info else None,
            'repeat_avg_len': repeat_info['repeat_avg_len'] if repeat_info else None,
            'spacer_avg_len': repeat_info['spacer_avg_len'] if repeat_info else None,
        }
        records.append(record)
    return records

def batch_process_crispr_files(df, crt_file_dir, crispr_distance=10000):
    updated_crispr_dfs = []

    for index, row in df.iterrows():
        target_file = row['target_file'].replace('.domtblout', '')
        source = row['source']

        # 生成crt文件对应的base_name
        if 'phage' in source:
            base_name = target_file.replace('.faa', '')
        elif '.fa.faa' in target_file:
            base_name = target_file.replace('.fa.faa', '')
        else:
            base_name = target_file.replace('.faa', '')

        crt_file_path = os.path.join(crt_file_dir, base_name + '.crt.txt')

        # 文件存在检查，避免出错
        if not os.path.isfile(crt_file_path):
            print(f"[警告] 缺失CRISPR文件：{crt_file_path}，跳过该条记录")
            continue

        crt_df = parse_crispr_file(crt_file_path)
        if len(crt_df) != 0:
            # genome字符串处理
            if 'IMGM_meta' in row['source']:
                genome = row['target_name'].values
            else:
                genome = "_".join(row['target_name'].split('_')[:-1])

            # 筛选子集，crt_df['genome']必须是genome的子串
            sub_df2 = crt_df[crt_df['genome'].apply(lambda x: x in genome)]
            sub_df1 = crt_df[crt_df['genome'] == genome]
            sub_df = pd.concat([sub_df1, sub_df2], ignore_index=True)
            sub_df = sub_df.applymap(lambda x: str(x) if isinstance(x, list) else x)
            sub_df = sub_df.drop_duplicates()

            crispr_coords = []
            print(len(sub_df))
            prot_start = row['prot_start']
            prot_end = row['prot_end']
            prot_id = row['target_name']
            window_start = prot_start - crispr_distance
            window_end = prot_end + crispr_distance

            for crispr_num, crispr_group in sub_df.groupby('crispr_num'):
                crispr_start = crispr_group['crispr_start'].min()
                crispr_end = crispr_group['crispr_end'].max()

                # 判断CRISPR是否在蛋白上下游crispr_distance范围内
                if window_start <= crispr_end and crispr_start <= window_end:
                    crispr_group = crispr_group.copy()
                    crispr_group['prot_start'] = prot_start
                    crispr_group['prot_end'] = prot_end
                    crispr_group['prot_id'] = prot_id

                    updated_crispr_dfs.append(crispr_group)
                    crispr_coords.append(f"{crispr_start}-{crispr_end}")

            if crispr_coords:
                crispr_summary = ";".join(crispr_coords)
                df.at[index, 'crispr_coords'] = crispr_summary
        else:
            df.at[index, 'crispr_coords'] = ''


    if updated_crispr_dfs:
        res_crispr_df = pd.concat(updated_crispr_dfs, ignore_index=True)
    else:
        res_crispr_df = pd.DataFrame()

    return df, res_crispr_df





if __name__ == '__main__':
    df = pd.read_table('/work/zhangrh/couple_procject/result/cas9/processed_df/step3_df.tsv')

    cluster_dir = "/work/zhangrh/couple_procject/result/cas9/cluster"  # 替换为聚类结果所在目录
    identity_thresholds = [100, 98, 90, 70, 50, 30]  # 可以根据实际需要修改
    identity_threshold = 50
    # 处理聚类结果并更新df
    df = process_clusters(cluster_dir, df, identity_thresholds)
    sub_df = df[df['target_name'] == df[f'cluster_{identity_threshold}_repseq']]
    sub_df.to_csv(f"/work/zhangrh/couple_procject/result/cas9/processed_df/step2_identity{identity_threshold}_g2.tsv", sep='\t', index=False)




    df.to_csv('/work/zhangrh/couple_procject/result/cas9/processed_df/step4_df.tsv', sep='\t', index=False)

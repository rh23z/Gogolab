import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed


def get_faa_seg_info(faa_file, target_seqid, source):
    if 'IMGM_metagenome' in source:
        faa_seg_info = process_gff_faa(faa_file, target_seqid)
    elif 'translated_cds' in faa_file and 'NCBI' in source:
        faa_seg_info = process_translated_CDS_faa(faa_file, target_seqid)
    else:
        faa_seg_info = process_prodigal_faa(faa_file, target_seqid)
    return faa_seg_info


def process_gff_faa(faa_path, target_seqid):
    gff_path = faa_path.replace('faa', 'gff')
    if not os.path.exists(gff_path):  # Check if GFF file exists
        print(f"Skipping {gff_path}, file does not exist.")
        return {}  # Return an empty dictionary if GFF file does not exist

    gff_df = pd.read_table(gff_path, header=None, comment='#')
    gff_df = gff_df.iloc[:, [0, 3, 4, 6, 8]]
    gff_df.columns = ['genome', 'start', 'end', 'strand', 'description']
    pattern = r"locus_tag=([^;]+)"
    gff_df['seqid'] = gff_df['description'].str.extract(pattern)

    target_genome = gff_df[gff_df['seqid'] == target_seqid]['genome'].unique()
    gff_df = gff_df[gff_df['genome'].isin(target_genome)]

    fasta_info = {}
    for record in SeqIO.parse(faa_path, "fasta"):
        rec_id = str(record.id)
        if rec_id in gff_df['seqid'].values:
            fasta_info[rec_id] = gff_df[gff_df['seqid'] == rec_id]['start'].values[0], \
            gff_df[gff_df['seqid'] == rec_id]['end'].values[0], \
                gff_df[gff_df['seqid'] == rec_id]['strand'].values[0], record
    return fasta_info


def process_prodigal_faa(faa_path, target_seqid):
    def extract_position_and_strand(record):
        description = record.description
        parts = description.split('#')
        if len(parts) < 4:
            raise ValueError(f"Description format error: {description}")
        start, end, strand = map(int, parts[1:4])
        return start, end, strand, record

    def get_fasta_info(fasta_file, target_seqid):
        fasta_info = {}
        target_genome = '_'.join(target_seqid.split('_')[:-1])
        for record in SeqIO.parse(fasta_file, "fasta"):
            rec_id = str(record.id)
            genome = '_'.join(rec_id.split('_')[:-1])  # 假设genome信息在ID的第一个部分

            if genome == target_genome:
                fasta_info[rec_id] = extract_position_and_strand(record)
        return fasta_info

    faa_seg_info = get_fasta_info(faa_path, target_seqid)
    return faa_seg_info


def process_translated_CDS_faa(faa_path, target_seqid):
    def extract_numbers(input_str):
        """
        提取字符串中的所有数字
        """
        numbers = re.findall(r'\d+', input_str)
        return list(map(int, numbers))

    def extract_position_and_strand2(record):
        description = record.description
        if 'complement' in description:
            strand = -1
        else:
            strand = 1
        parts = description.split('[location=')[-1].split(']')[0].split(',')[0]
        positions = extract_numbers(parts)

        sorted_positions = np.sort(positions)  # 对数字进行排序
        start = int(sorted_positions[0])
        end = int(sorted_positions[-1])

        return start, end, strand, record

    def get_fasta_info(fasta_file, target_seqid):
        fasta_info = {}
        target_genome = target_seqid.replace('lcl|', '').split('_prot_')[0]
        for record in SeqIO.parse(fasta_file, "fasta"):
            rec_id = str(record.id)
            genome = record.id.replace('lcl|', '').split('_prot_')[0]
            if genome == target_genome:
                fasta_info[rec_id] = extract_position_and_strand2(record)

        return fasta_info

    faa_seg_info = get_fasta_info(faa_path, target_seqid)
    return faa_seg_info


# 根据目标蛋白上下游的起止位置，提取对应片段信息
def extract_flanking_segments(faa_seg_info, target_name, target_start, target_end, upstream=10000, downstream=10000):
    flanking_segments = []
    flanking_segments_info = []
    # 计算上下游范围
    upstream_start = target_start - upstream
    downstream_end = target_end + downstream

    for record_id, segment in faa_seg_info.items():
        segment_start, segment_end, segment_strand, record = segment
        # 判断该片段是否在上下游范围内
        if (segment_end >= upstream_start and segment_start <= downstream_end):
            # 记录该片段的FASTA ID和起止位置
            flanking_segments_info.append(
                (record.id, segment_start, segment_end, segment_strand)
            )
            flanking_segments.append(record)
    return flanking_segments_info, flanking_segments


# 提取片段的FASTA序列并写入到新的FASTA文件
def write_flanking_segments_to_fasta(flanking_segments, output_fasta):
    # 创建新的FASTA文件
    with open(output_fasta, 'w') as fasta_file:
        for record in flanking_segments:
            # 创建一个新的 SeqRecord 对象，并写入文件
            fasta_record = SeqRecord(record.seq, id=record.id, description=record.description)
            SeqIO.write(fasta_record, fasta_file, "fasta")


# 读取DataFrame并处理
def process_single_row(index, row, faa_dict):
    fasta_file_name = row['target_file'].split('.domtblout')[0]
    root_faa_dir = faa_dict[row['source']]
    faa_file = os.path.join(root_faa_dir, fasta_file_name)

    id_ = row['target_name']
    target_start, target_end, target_strand, source = min(row['prot_start'], row['prot_end']), max(
        row['prot_start'], row['prot_end']), row['strand'], row['source']

    faa_seg_info = get_faa_seg_info(faa_file, id_, source)

    # If faa_seg_info is empty, return empty lists
    if not faa_seg_info:
        return index, [], []

    flanking_segments_info, flanking_segments = extract_flanking_segments(faa_seg_info, id_, target_start, target_end)

    # 返回处理后的行数据
    return index, flanking_segments_info, flanking_segments


def process_dataframe_parallel(df, out_dir, faa_dict):
    all_flanking_segments = {}
    futures = []

    with ProcessPoolExecutor() as executor:
        for index, row in df.iterrows():
            future = executor.submit(process_single_row, index, row, faa_dict)
            futures.append(future)

        for future in as_completed(futures):
            index, flanking_segments_info, flanking_segments = future.result()

            # 将片段信息记录到新列
            df.loc[index, 'flanking_segments'] = str(flanking_segments_info)

            # 将片段添加到集合中
            all_flanking_segments[index] = flanking_segments

    # 将所有片段写入FASTA文件
    for index, segments in all_flanking_segments.items():
        output_fasta = os.path.join(out_dir, f"{df.loc[index, 'target_name']}_flanking.fasta")
        write_flanking_segments_to_fasta(segments, output_fasta)

    return df


if __name__ == '__main__':
    ids = [
        'archaea_assembly_summary',
        'bacteria_assembly_summary',
        'BGItrans',
        'metagenomes_assembly_summary',
        'phagescope',
        'viral_assembly_summary',
        'CRBC',
        'IMGM',
        'MGnify',
        'IMGM_metagenome'
    ]

    faa_dict = {
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

    # 读取DataFrame并处理
    df = pd.read_table(
        '/work/zhangrh/couple_procject/result/cas9/raw_data/step1_df_all_function_filtered.tsv')  # 替换为实际DataFrame路径

    # 处理DataFrame
    tags = ids
    filtered_df = df[df['source'].apply(lambda x: any(tag in x for tag in tags))]
    final_filtered_df = process_dataframe_parallel(filtered_df, '/work/zhangrh/couple_procject/result/cas9/fasta/step3',
                                                   faa_dict)

    # 输出新的DataFrame
    final_filtered_df.to_csv('/work/zhangrh/couple_procject/result/cas9/processed_df/step3_df.tsv', sep='\t',
                             index=False)

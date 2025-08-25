import os
import pandas as pd
from Bio import SeqIO
import gzip
from concurrent.futures import ProcessPoolExecutor, as_completed

# 配置
output_dir = '/work/zhangrh/couple_procject/result/cas9/fasta/step4_sup'
log_file_path = '/work/zhangrh/couple_procject/result/cas9/array/sup.log'
check_dir = '/work/zhangrh/couple_procject/result/cas9/array/crt'  # 用于检查已有文件的文件夹

source_dict = {
    'archaea_assembly_summary': "/work/data_share/NCBI_zrh/archaea_assembly_summary/",
    'bacteria_assembly_summary': "/work/data_share/NCBI_zrh/bacteria_assembly_summary/",
    'BGItrans': "/work/data/BGItrans/final_bins/",
    'metagenomes_assembly_summary': "/work/data_share/NCBI_zrh/metagenomes_assembly_summary/",
    'phagescope': '/work/data/phage/phagescope_all/',
    'viral_assembly_summary': "/work/data_share/NCBI_zrh/viral_assembly_summary/",
    'IMGM': '/work/data/IMGM_1025/fna/',
    'MGnify': '/work/data/MGnify/fna/',
    'CRBC': '/work/data/CRBC_genome/fna/',
    'IMGM_metagenome': '/work/data/IMGM_metagenome/fna/'
}

# 读取表格
step1_df = pd.read_table(r'/work/zhangrh/couple_procject/result/cas9/raw_data/step1_df_all_function_filtered_g2.tsv')
step1_df = step1_df[step1_df['source']=='IMGM_metagenome']
def decompress_gz_to_tmp(gz_file, suffix=".fna"):
    """解压.gz文件到临时文件并返回路径"""
    tmp_file = f"{gz_file}_tmp{suffix}"
    with gzip.open(gz_file, 'rb') as f_in:
        with open(tmp_file, 'wb') as f_out:
            f_out.write(f_in.read())
    return tmp_file
import re
from threading import Lock

# ========== GFF 缓存 ==========
_GFF_CACHE = {}
_GFF_LOCK = Lock()

def load_gff_two_cols(gff_path: str) -> pd.DataFrame:
    """只读需要的两列并缓存：genome(第0列), description(第8列)"""
    with _GFF_LOCK:
        df_cached = _GFF_CACHE.get(gff_path)
        if df_cached is not None:
            return df_cached

    df = pd.read_table(
        gff_path,
        sep='\t',
        header=None,
        comment='#',
        usecols=[0, 8],
        names=['genome', 'description'],
        dtype={'genome': 'string', 'description': 'string'},
        engine='c',
        memory_map=True
    )
    # 压缩一下 genome 列，降低常驻内存
    df['genome'] = df['genome'].astype('category')

    with _GFF_LOCK:
        _GFF_CACHE[gff_path] = df
    return df

def process_target_file(group_df):
    """处理同一个 target_file 的所有行，输出一个文件（加速 GFF 处理版）"""
    row0 = group_df.iloc[0]
    target_file = row0['target_file'].replace('.domtblout', '')
    source = row0['source']
    source_dir = source_dict[source]

    if 'phage' in source_dir:
        base_name = target_file.replace('.faa', '')
        target_file_name = base_name + '.fasta'
    elif '.fa.faa' in target_file:
        base_name = target_file.replace('.fa.faa', '')
        target_file_name = base_name + '.fa'
    else:
        base_name = target_file.replace('.faa', '')
        target_file_name = base_name + '.fna'

    # 先检查 check_dir 是否有文件包含 base_name
    try:
        if os.path.exists(check_dir,base_name+'.crt.txt'):
            return f"[跳过] 文件夹 {check_dir} 已存在包含 {base_name} 的文件"
    except Exception:
        # check_dir 不存在或无法读取时，按无匹配处理
        pass

    fna_path_unz = os.path.join(source_dir, target_file_name)
    fna_path_gz = fna_path_unz + ".gz"
    output_file = os.path.join(output_dir, target_file_name)

    if os.path.exists(output_file):
        return f"[跳过] 输出文件已存在: {output_file}"

    # ===================== 收集 genome_ids（加速） =====================
    genome_ids_set = set()

    # IMGM_metagenome 通过 GFF 解析
    imgm_mask = group_df['source'].str.contains('IMGM_meta', na=False)
    if imgm_mask.any():
        gff_path = fna_path_unz.replace('fna', 'gff')
        try:
            gff_df = load_gff_two_cols(gff_path)  # 只读 0/8 两列并缓存
        except Exception as e:
            msg = f"[失败] 读取GFF失败: {gff_path}，错误：{e}\n"
            with open(log_file_path, 'a') as log_f:
                log_f.write(msg)
            return f"[失败] {target_file_name} 读取GFF失败，错误已记录"

        # 本组内唯一 target_name
        unique_targets = group_df.loc[imgm_mask, 'target_name'].dropna().unique().tolist()

        # 分块避免正则过长
        CHUNK = 200
        for i in range(0, len(unique_targets), CHUNK):
            chunk = unique_targets[i:i+CHUNK]
            try:
                pat = '|'.join(re.escape(t) for t in chunk)
                candidate = gff_df[gff_df['description'].str.contains(pat, regex=True, na=False)]
            except re.error:
                candidate = gff_df  # 退化为全表候选

            # 对每个 target 再做精确 contains（regex=False 避免特殊字符误伤）
            for t in chunk:
                matched = candidate.loc[
                    candidate['description'].str.contains(t, regex=False, na=False),
                    'genome'
                ]
                if not matched.empty:
                    for g in matched.astype(str).unique().tolist():
                        if g:
                            genome_ids_set.add(g)

    # 其它来源：直接从 target_name 推断
    for _, row in group_df.loc[~imgm_mask].iterrows():
        genome = "_".join(row['target_name'].split('_')[:-1])
        if genome:
            genome_ids_set.add(genome)

    print(f"[处理] {target_file_name}，包含 {len(genome_ids_set)} 个序列候选")

    # ===================== 打开输入文件（原逻辑不变） =====================
    if os.path.exists(fna_path_unz):
        input_path = fna_path_unz
    elif os.path.exists(fna_path_gz):
        try:
            input_path = decompress_gz_to_tmp(fna_path_gz)
        except Exception as e:
            msg = f"[失败] 解压失败: {fna_path_gz}，错误：{e}\n"
            with open(log_file_path, 'a') as log_f:
                log_f.write(msg)
            return f"[失败] {target_file_name} 解压失败，错误已记录"
    else:
        msg = f"[未找到] {fna_path_unz} 或 {fna_path_gz}\n"
        with open(log_file_path, 'a') as log_f:
            log_f.write(msg)
        return msg.strip()

    # ===================== 扫描 FASTA 并写出（仅改为 set 查找） =====================
    try:
        records_to_write = []
        with open(input_path, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                if record.id in genome_ids_set:
                    records_to_write.append(record)

        if records_to_write:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            with open(output_file, 'w') as out_f:
                SeqIO.write(records_to_write, out_f, 'fasta')
            print(f"[成功] 写入 {len(records_to_write)} 条序列到 {output_file}")
            return f"[成功] 写入 {len(records_to_write)} 条序列到 {output_file}"
        else:
            msg = f"[警告] 没有匹配的序列在 {target_file_name} 中"
            with open(log_file_path, 'a') as log_f:
                log_f.write(msg + "\n")
            return msg
    except Exception as e:
        msg = f"[失败] {target_file_name} 读取或写入失败，错误：{e}\n"
        with open(log_file_path, 'a') as log_f:
            log_f.write(msg)
        return f"[失败] {target_file_name} 处理失败，错误已记录"

# 按 target_file 分组并处理
from concurrent.futures import ThreadPoolExecutor, as_completed
import traceback

def parallel_process(df, max_workers=32):
    """
    使用线程池并行处理 step1_df 按 target_file 分组的任务（带 tqdm 进度条）
    """
    from tqdm.auto import tqdm  # 放在函数内部，避免全局改动
    results = []

    # 按 target_file 分组，生成独立副本
    grouped = [(target_file, group_df.copy()) for target_file, group_df in df.groupby('target_file')]
    total = len(grouped)

    if total == 0:
        return results

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # 提交所有任务
        future_to_group = {executor.submit(process_target_file, group): target_file for target_file, group in grouped}

        # 进度条
        with tqdm(total=total, desc="处理 target_file 组", unit="file", dynamic_ncols=True) as pbar:
            for future in as_completed(future_to_group):
                target_file = future_to_group[future]
                try:
                    result = future.result()
                except Exception as e:
                    tb = traceback.format_exc()
                    result = f"[失败] {target_file} 线程处理失败，错误：{e}\n{tb}"

                # 进度安全输出 & 更新进度
                pbar.write(result)
                results.append(result)
                pbar.update(1)
                # 在进度条尾部展示最近处理的 target_file（截断显示更清晰）
                pbar.set_postfix_str(str(target_file)[-40:])

    return results



# 执行并行处理
parallel_results = parallel_process(step1_df, max_workers=16)  # 可根据机器调整线程数
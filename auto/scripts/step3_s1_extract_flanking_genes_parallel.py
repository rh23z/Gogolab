from __future__ import annotations
import os
import re
import sys
import json
import gzip
import argparse
import numpy as np
import pandas as pd
from typing import Dict, Tuple, Iterable, List
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ------------------------------
# Helpers for IO & parsing
# ------------------------------

def _smart_open_fasta(path: str):
    """Return a text handle for FASTA path (supports .gz)."""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def _iter_seqrecords(path: str) -> Iterable[SeqRecord]:
    with _smart_open_fasta(path) as handle:
        for rec in SeqIO.parse(handle, 'fasta'):
            yield rec


def _possible_faa_paths(root_dir: str, base: str) -> List[str]:
    """Return candidate FASTA paths given a root dir and basename without extension.
    We try common extensions and their .gz variants.
    """
    exts = ['.faa', '.fasta', '.fa']
    candidates = []
    for ext in exts:
        p = os.path.join(root_dir, base + ext)
        candidates.append(p)
        candidates.append(p + '.gz')
    return candidates


def _resolve_faa_path(root_dir: str, domtblout_basename: str) -> str | None:
    """domtblout_basename is the filename without trailing ".domtblout".
    We search for <basename>.faa/.fasta/.fa(.gz) under root_dir.
    """
    base = domtblout_basename
    for cand in _possible_faa_paths(root_dir, base):
        if os.path.exists(cand):
            return cand
    return None


# ------------------------------
# Source-specific parsers
# ------------------------------

def process_gff_faa(faa_path: str, target_seqid: str) -> Dict[str, Tuple[int, int, int, SeqRecord]]:
    """For IMGM_metagenome-like sources with paired GFF.
    We infer the GFF path by swapping extension.
    """
    # Derive gff path robustly
    if faa_path.endswith('.gz'):
        gff_path = re.sub(r'\.(fa|faa|fasta)\.gz$', '.gff', faa_path)
    else:
        gff_path = re.sub(r'\.(fa|faa|fasta)$', '.gff', faa_path)

    if not os.path.exists(gff_path):
        raise FileNotFoundError(f"GFF not found for {faa_path}: {gff_path}")

    gff_df = pd.read_table(gff_path, header=None, comment='#')
    gff_df = gff_df.iloc[:, [0, 3, 4, 6, 8]].copy()
    gff_df.columns = ['genome', 'start', 'end', 'strand', 'description']
    pattern = r"locus_tag=([^;]+)"
    gff_df['seqid'] = gff_df['description'].str.extract(pattern)

    target_genome = gff_df.loc[gff_df['seqid'] == target_seqid, 'genome'].unique()
    gff_df = gff_df[gff_df['genome'].isin(target_genome)]

    fasta_info: Dict[str, Tuple[int, int, int, SeqRecord]] = {}
    for record in _iter_seqrecords(faa_path):
        rec_id = str(record.id)
        if rec_id in set(gff_df['seqid']):
            sub = gff_df[gff_df['seqid'] == rec_id]
            # 预期只有一行；若多行，取第一行
            start = int(sub['start'].iloc[0])
            end = int(sub['end'].iloc[0])
            strand = int(1 if str(sub['strand'].iloc[0]) == '+' else -1)
            fasta_info[rec_id] = (start, end, strand, record)
    return fasta_info


def process_prodigal_faa(faa_path: str, target_seqid: str) -> Dict[str, Tuple[int, int, int, SeqRecord]]:
    def extract_position_and_strand(record: SeqRecord):
        description = record.description
        parts = description.split('#')
        if len(parts) < 4:
            raise ValueError(f"Description format error: {description}")
        start, end, strand = map(int, parts[1:4])
        return start, end, strand, record

    target_genome = '_'.join(target_seqid.split('_')[:-1])
    fasta_info: Dict[str, Tuple[int, int, int, SeqRecord]] = {}

    for record in _iter_seqrecords(faa_path):
        rec_id = str(record.id)
        genome = '_'.join(rec_id.split('_')[:-1])
        if genome == target_genome:
            fasta_info[rec_id] = extract_position_and_strand(record)
    return fasta_info


def process_translated_CDS_faa(faa_path: str, target_seqid: str) -> Dict[str, Tuple[int, int, int, SeqRecord]]:
    def extract_numbers(input_str: str) -> List[int]:
        numbers = re.findall(r'\d+', input_str)
        return list(map(int, numbers))

    def extract_position_and_strand2(record: SeqRecord):
        description = record.description
        strand = -1 if 'complement' in description else 1
        parts = description.split('[location=')[-1].split(']')[0].split(',')[0]
        positions = extract_numbers(parts)
        sorted_positions = np.sort(positions)
        start = int(sorted_positions[0])
        end = int(sorted_positions[-1])
        return start, end, strand, record

    target_genome = target_seqid.replace('lcl|', '').split('_prot_')[0]
    fasta_info: Dict[str, Tuple[int, int, int, SeqRecord]] = {}

    for record in _iter_seqrecords(faa_path):
        genome = record.id.replace('lcl|', '').split('_prot_')[0]
        if genome == target_genome:
            fasta_info[str(record.id)] = extract_position_and_strand2(record)
    return fasta_info


def get_faa_seg_info(faa_file: str, target_seqid: str, source: str) -> Dict[str, Tuple[int, int, int, SeqRecord]]:
    if 'IMGM_metagenome(no more use)' in source:
        return process_gff_faa(faa_file, target_seqid)
    elif 'translated_cds' in faa_file and 'NCBI' in source:
        return process_translated_CDS_faa(faa_file, target_seqid)
    else:
        return process_prodigal_faa(faa_file, target_seqid)


# ------------------------------
# Flanking extraction
# ------------------------------

def extract_flanking_segments(
    faa_seg_info: Dict[str, Tuple[int, int, int, SeqRecord]],
    target_name: str,
    target_start: int,
    target_end: int,
    upstream: int = 10000,
    downstream: int = 10000,
):
    flanking_segments: List[SeqRecord] = []
    flanking_segments_info: List[Tuple[str, int, int, int]] = []

    upstream_start = target_start - upstream
    downstream_end = target_end + downstream

    for record_id, segment in faa_seg_info.items():
        segment_start, segment_end, segment_strand, record = segment
        if (segment_end >= upstream_start and segment_start <= downstream_end):
            flanking_segments_info.append((record.id, int(segment_start), int(segment_end), int(segment_strand)))
            flanking_segments.append(record)
    return flanking_segments_info, flanking_segments


def write_flanking_segments_to_fasta(flanking_segments: List[SeqRecord], output_fasta: str):
    if not flanking_segments:
        # 没有邻近片段也要产生一个空文件以便后续检查（可选）
        open(output_fasta, 'w').close()
        return
    with open(output_fasta, 'w') as fasta_file:
        for record in flanking_segments:
            fasta_record = SeqRecord(record.seq, id=record.id, description=record.description)
            SeqIO.write(fasta_record, fasta_file, 'fasta')


# ------------------------------
# Row processing (for parallel)
# ------------------------------

def _process_single_row(index: int, row: pd.Series, faa_roots: Dict[str, str], upstream: int, downstream: int,
                         overwrite: bool, out_dir: str):
    # domtblout 基名 -> 原始 FASTA 基名
    fasta_base = str(row['target_file']).split('.domtblout')[0]
    source = str(row['source'])
    root_faa_dir = faa_roots.get(source)
    if root_faa_dir is None:
        return index, '[]', []  # unknown source -> skip

    faa_file = _resolve_faa_path(root_faa_dir, fasta_base)
    if faa_file is None:
        return index, '[]', []  # missing fasta -> skip

    id_ = row['target_name']
    target_start = int(min(row['prot_start'], row['prot_end']))
    target_end = int(max(row['prot_start'], row['prot_end']))

    try:
        faa_seg_info = get_faa_seg_info(faa_file, id_, source)
        flanking_segments_info, flanking_segments = extract_flanking_segments(
            faa_seg_info, id_, target_start, target_end, upstream=upstream, downstream=downstream
        )
    except Exception as e:
        # 返回空并在主进程里打印
        return index, f"ERROR:{e}", []

    # 写 FASTA（每个 target_name 一份）
    out_fa = os.path.join(out_dir, f"{id_}_flanking.fasta")
    if overwrite or (not os.path.exists(out_fa)):
        os.makedirs(out_dir, exist_ok=True)
        write_flanking_segments_to_fasta(flanking_segments, out_fa)

    return index, json.dumps(flanking_segments_info), []  # sequences 已写文件，这里不回传


def process_dataframe_parallel(df: pd.DataFrame, out_dir: str, faa_roots: Dict[str, str],
                               workers: int, upstream: int, downstream: int, overwrite: bool) -> pd.DataFrame:
    df = df.copy()
    df['flanking_segments'] = ''

    futures = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        for index, row in df.iterrows():
            futures.append(
                executor.submit(_process_single_row, index, row, faa_roots, upstream, downstream, overwrite, out_dir)
            )
        for i, fut in enumerate(as_completed(futures), 1):
            index, info_json, _ = fut.result()
            if info_json.startswith('ERROR:'):
                print(f"[WARN] index={index} {info_json}")
                df.loc[index, 'flanking_segments'] = '[]'
            else:
                df.loc[index, 'flanking_segments'] = info_json
            if i % 100 == 0:
                print(f"… processed {i}/{len(futures)} rows")

    return df


# ------------------------------
# CLI
# ------------------------------

def _parse_source_map(args) -> Dict[str, str]:
    # Priority: --source-map file > repeated --map k=v > built-in defaults
    if args.source_map:
        path = args.source_map
        try:
            # Try JSON first
            with open(path, 'r') as f:
                txt = f.read().strip()
                try:
                    return json.loads(txt)
                except json.JSONDecodeError:
                    pass
            # Try YAML (optional)
            try:
                import yaml  # type: ignore
                with open(path, 'r') as f:
                    return yaml.safe_load(f)
            except Exception:
                raise ValueError(f"无法解析 --source-map: {path} (支持 JSON/YAML)")
        except FileNotFoundError:
            raise

    if args.map:
        mapping: Dict[str, str] = {}
        for item in args.map:
            if '=' not in item:
                raise ValueError(f"--map 需要 key=val 形式，得到: {item}")
            k, v = item.split('=', 1)
            mapping[k.strip()] = v.strip()
        return mapping

    # built-in defaults (与历史版本一致)
    return {
        'archaea_assembly_summary': "/work/data_share/NCBI_zrh/archaea_assembly_summary_prodigal/",
        'bacteria_assembly_summary': "/work/data_share/NCBI_zrh/bacteria_assembly_summary_prodigal/",
        'BGItrans': "/work/data/BGItrans/final_bins_faa/",
        'metagenomes_assembly_summary': "/work/data_share/NCBI_zrh/metagenomes_assembly_summary_prodigal/",
        'phagescope': '/work/data/phage/phagescope_all_prodigal/',
        'viral_assembly_summary': "/work/data_share/NCBI_zrh/viral_assembly_summary_prodigal/",
        'IMGM': "/work/data/IMGM_metagenome/faa/",
        'CRBC': "/work/data/CRBC_genome/faa/",
        'IMGM_metagenome': "/work/data/IMGM_metagenome/prodigal/",
        'MGnify': "/work/data/MGnify/faa/",

    }


def main():
    p = argparse.ArgumentParser(
        description='Extract flanking neighbour proteins in parallel (Step 3).'
    )
    p.add_argument('-i', '--summary', required=True, help='输入 summary.tsv (包含 target_name/target_file/prot_start/prot_end/strand/source 等列)')
    p.add_argument('-o', '--out-dir', required=True, help='输出 FASTA 的目录（每个 target_name 一份）')
    p.add_argument('--df-out', default=None, help='处理后的 DataFrame 输出路径（默认 <out-dir>/step3_df.tsv）')
    p.add_argument('--source-map', default=None, help='source 映射 JSON/YAML 文件路径；或使用 --map 重复指定')
    p.add_argument('--map', action='append', default=None, help='单个 source 映射，形如 key=/abs/path （可重复）')
    p.add_argument('--filter-sources', default=None, help='仅处理这些 source（逗号分隔）；默认处理映射中所有 key')
    p.add_argument('--workers', type=int, default=os.cpu_count() or 8, help='并行进程数（默认=CPU核数）')
    p.add_argument('--upstream', type=int, default=10000, help='上游窗口大小(bp)')
    p.add_argument('--downstream', type=int, default=10000, help='下游窗口大小(bp)')
    p.add_argument('--overwrite', action='store_true', help='覆盖已存在的 <target_name>_flanking.fasta')
    args = p.parse_args()

    # Load table
    df = pd.read_table(args.summary)

    # Build mapping
    faa_roots = _parse_source_map(args)

    # Filter by sources if requested
    if args.filter_sources:
        keep = set([s.strip() for s in args.filter_sources.split(',') if s.strip()])
    else:
        keep = set(faa_roots.keys())



    # Ensure out-dir exists
    os.makedirs(args.out_dir, exist_ok=True)

    # Process
    print(f"[INFO] rows={len(df)} workers={args.workers} upstream={args.upstream} downstream={args.downstream}")
    df2 = process_dataframe_parallel(
        df,
        out_dir=args.out_dir,
        faa_roots=faa_roots,
        workers=args.workers,
        upstream=args.upstream,
        downstream=args.downstream,
        overwrite=args.overwrite,
    )

    # Save DF
    df_out = args.df_out or os.path.join(args.out_dir, 'step3_df.tsv')
    os.makedirs(os.path.dirname(df_out), exist_ok=True)
    df2.to_csv(df_out, sep='\t', index=False)
    print(f"[OK] DataFrame saved -> {df_out}")


if __name__ == '__main__':
    main()

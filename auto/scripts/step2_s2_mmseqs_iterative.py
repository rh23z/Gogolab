#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path
from datetime import datetime

def log(msg):
    print(f"[{datetime.now().strftime('%F %T')}] {msg}")

def run_mmseqs(input_fasta, output_prefix, tmp_dir, min_seq_id, cov, threads):
    cmd = [
        "mmseqs", "easy-cluster", str(input_fasta), str(output_prefix),
        str(tmp_dir),
        "--cov-mode", "0",
        "-c", str(cov),
        "--min-seq-id", f"{min_seq_id:.2f}",
        "--cluster-reassign", "1",
        "--threads", str(threads)
    ]
    log(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def find_rep_seq_file(prefix):
    rep_fasta = Path(f"{prefix}_rep_seq.fasta")
    rep_faa = Path(f"{prefix}_rep_seq.faa")
    if rep_fasta.exists():
        return rep_fasta
    elif rep_faa.exists():
        return rep_faa
    else:
        raise FileNotFoundError(f"Cannot find representative sequence file for prefix {prefix}")

def main():
    parser = argparse.ArgumentParser(description="Iterative MMseqs clustering in Python")
    parser.add_argument("input_fasta", help="Initial input FASTA file")
    parser.add_argument("output_dir", help="Output directory to save clustering results")
    parser.add_argument("tmp_dir", help="Temporary directory for MMseqs")
    parser.add_argument("--identity", nargs="+", type=int, default=[100, 98, 90, 70, 50, 30],
                        help="List of identity thresholds (in %), e.g. 100 98 90")
    parser.add_argument("-c", "--coverage", type=float, default=0.8,
                        help="Coverage threshold (default: 0.8)")
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="Number of threads for MMseqs")

    args = parser.parse_args()
    input_fasta = Path(args.input_fasta).resolve()
    output_dir = Path(args.output_dir).resolve()
    tmp_dir = Path(args.tmp_dir).resolve()

    output_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    current_input = input_fasta
    for round_num, id_pct in enumerate(args.identity, start=1):
        min_id = id_pct / 100
        prefix = output_dir / f"round{round_num}_{id_pct}"
        log(f"==> Round {round_num}: identity={min_id:.2f}, coverage={args.coverage}")
        run_mmseqs(current_input, prefix, tmp_dir, min_id, args.coverage, args.threads)
        current_input = find_rep_seq_file(prefix)

    log("✅ All clustering rounds completed.")

if __name__ == "__main__":
    main()

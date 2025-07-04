# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import os
import logging
from datetime import datetime
import pandas as pd
from Bio import SeqIO

def run_metrics(args):

    def calculate_n50_l50(lengths):
        lengths_sorted = sorted(lengths, reverse=True)
        total = sum(lengths_sorted)
        cumulative = 0
        for i, length in enumerate(lengths_sorted):
            cumulative += length
            if cumulative >= total / 2:
                return length, i + 1
        return 0, 0
    
    def compute_fasta_metrics(filepath):
        records = list(SeqIO.parse(filepath, "fasta"))
        if not records:
            return None
        
        lengths = [len(record.seq) for record in records]
        total_length = sum(lengths)
        n_count = sum(str(record.seq).upper().count('N') for record in records)
        gc_count = sum(str(record.seq).upper().count('G') + str(record.seq).upper().count('C') for record in records)
        gc_content = gc_count / total_length * 100 if total_length > 0 else 0
        n50, l50 = calculate_n50_l50(lengths)

        return {
            "Filename": os.path.basename(filepath),
            "Num_Contigs": len(lengths),
            "Total_Length": total_length,
            "Longest_Contig": max(lengths),
            "Shortest_Contig": min(lengths),
            "GC_Content(%)": round(gc_content, 2),
            "N_Count": n_count,
            "N50": n50,
            "L50": l50,
            "Contig_Lengths": ','.join(map(str, sorted(lengths, reverse=True)))
        }
    
    logger = logging.getLogger("metrics")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if not any(isinstance(h, logging.FileHandler) and h.baseFilename.endswith("seqforge_metrics.log") for h in logger.handlers):
        file_handler = logging.FileHandler("seqforge_metrics.log", mode='a')
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    start_time = datetime.now()
    print(f"Metrics analysis started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Metrics analysis started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    input_path = os.path.abspath(args.fasta_directory)
    if not os.path.exists(input_path):
        msg = f"Error: The specified input does not exist: {args.fasta_directory}"
        print(f"\n\033[91m{msg}\033[0m")
        logger.warning(msg)
        return
    
    fasta_files = []
    if os.path.isfile(input_path):
        fasta_files.append(input_path)
    else:
        for f in os.listdir(input_path):
            if f.endswith((".fasta", ".fa", ".fna", ".ffn", ".fas", ".faa")):
                fasta_files.append(os.path.join(input_path, f))

    if not fasta_files:
        msg = "No FASTA files found for processing."
        print(f"\n\033[93m{msg}\033[0m")
        logger.warning(msg)
        return
    
    metrics = []
    for file in fasta_files:
        result = compute_fasta_metrics(file)
        if result:
            metrics.append(result)

    df = pd.DataFrame(metrics)
    print_df = df.drop(columns=['Contig_Lengths'], errors='ignore')
    print("\n\033[92mFASTA Metrics Summary:\033[0m")
    print(print_df.to_string(index=False))

    output_file = args.output or "fasta_metrics_summary.csv"
    df.to_csv(output_file, index=False)
    print(f"\n\033[92mSummary written to: {output_file}\033[0m")
    logger.info(f"Summary written to: {output_file}")

    end_time = datetime.now()
    logger.info(f"Metrics analysis completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    msg = f"Total runtime: {str(end_time - start_time)}"
    logger.info(msg)
    print(msg)

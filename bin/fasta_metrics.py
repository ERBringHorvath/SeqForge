# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import os
import logging
from datetime import datetime
import pandas as pd
from Bio import SeqIO

def run_metrics(args):

    def calculate_nx_lx(lengths, threshold=0.5):
        lengths_sorted = sorted(lengths, reverse=True)
        total = sum(lengths_sorted)
        cumulative = 0
        for i, length in enumerate(lengths_sorted):
            cumulative += length
            if cumulative >= total * threshold:
                return length, i + 1
        return 0, 0
    
    def compute_fasta_metrics(filepath, min_contig_size=500):
        records = list(SeqIO.parse(filepath, "fasta"))
        if not records:
            return None
        
        lengths = [len(record.seq) for record in records if len(record.seq) >= min_contig_size]
        if not lengths:
            return {
                "Filename": os.path.basename(filepath),
                "Num_Contigs": 0,
                "Num_Contigs_≥1kb": 0,
                "Num_Contigs_≥10kb": 0,
                "Num_Contigs_≥50kb": 0,
                "Num_Contigs_≥100kb": 0,
                "Total_Length": 0,
                "Total_Length_≥1kb": 0,
                "Longest_Contig": 0,
                "Shortest_Contig": 0,
                "GC_Content(%)": 0,
                "N_Count": 0,
                "N50": 0,
                "L50": 0,
                "N90": 0,
                "L90": 0,
                "Contig_Lengths": ''
            }
        
        total_length = sum(lengths)
        n_count = sum(str(record.seq).upper().count('N') for record in records if len(record.seq) >= min_contig_size)
        gc_count = sum(str(record.seq).upper().count('G') + str(record.seq).upper().count('C')
                       for record in records if len(record.seq) >= min_contig_size)
        gc_content = gc_count / total_length * 100 if total_length > 0 else 0

        n50, l50 = calculate_nx_lx(lengths, threshold=0.5)
        n90, l90 = calculate_nx_lx(lengths, threshold=0.9)

        contigs_1kb = sum(1 for l in lengths if l >= 1000)
        contigs_10kb = sum(1 for l in lengths if l >= 10000)
        contigs_50kb = sum(1 for l in lengths if l >= 50000)
        contigs_100kb = sum(1 for l in lengths if l >= 100000)
        total_length_1kb = sum(l for l in lengths if l >= 1000)

        contig_lengths_str = '|'.join(str(length) for length in sorted(lengths, reverse=True))

        return {
            "Filename": os.path.basename(filepath),
            "Num_Contigs": len(lengths),
            "Num_Contigs_≥1kb": contigs_1kb,
            "Num_Contigs_≥10kb": contigs_10kb,
            "Num_Contigs_≥50kb": contigs_50kb,
            "Num_Contigs_≥100kb": contigs_100kb,
            "Total_Length": total_length,
            "Total_Length_≥1kb": total_length_1kb,
            "Longest_Contig": max(lengths),
            "Shortest_Contig": min(lengths),
            "GC_Content(%)": round(gc_content, 2),
            "N_Count": n_count,
            "N50": n50,
            "N90": n90,
            "L50": l50,
            "L90": l90,
            "Contig_Lengths": contig_lengths_str
        }

    logger = logging.getLogger("metrics")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if not any(isinstance(h, logging.FileHandler) and h.baseFilename.endswith("seqforge_metrics.log") for h in logger.handlers):
        file_handler = logging.FileHandler("seqforge_metrics.log", mode='a')
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    logger.info(f"Metrics calculated based on minimum contig size of {args.min_contig_size} bp")
    print(f"\n\033[96mMetrics calculated based on minimum contig size of {args.min_contig_size} bp\033[0m")

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
        result = compute_fasta_metrics(file, min_contig_size=args.min_contig_size)
        if result:
            metrics.append(result)

    df = pd.DataFrame(metrics)
    print_df = df.drop(columns=['Contig_Lengths', 'Num_Contigs_≥1kb', 'Num_Contigs_≥10kb', 'Num_Contigs_≥50kb', 'Num_Contigs_≥100kb'], errors='ignore')
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
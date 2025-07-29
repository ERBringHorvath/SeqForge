# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import os
from datetime import datetime
import logging
import pandas as pd
from Bio import SeqIO
import concurrent.futures
from threading import Lock

from utils.file_handler import collect_fasta_files, cleanup_temp_dir

def process_contig_entry(row, fasta_map, evalue, min_perc, min_cov,
                         extracted_contigs, lock, *, logger):
    
    if float(row['evalue']) > evalue or float(row['pident']) < min_perc or float(row['query_coverage']) < min_cov:
        msg = (f"Skipping sequence ({row['database']} query: {row['query_file_name']} "
               f"pident: {row['pident']} query coverage: {row['query_coverage']} "
               f"evalue: {row['evalue']}) due to filtering thresholds")
        print(f"\033[91m{msg}\033[0m")
        logger.info(msg)
        return None

    contig_key = (row['database'], row['sseqid'])

    #Skip duplicates
    with lock:
        if contig_key in extracted_contigs:
            return None
        extracted_contigs.add(contig_key)

    original_fasta = fasta_map.get(row['database'])
    if not original_fasta:
        msg = f"No matching FASTA file found for {row['database']}"
        print(f"\033[91m{msg}\033[0m")
        logger.warning(msg)
        return None

    found = False
    for seq_record in SeqIO.parse(original_fasta, "fasta"):
        if seq_record.id == row['sseqid']:
            found = True
            header_id = f"{seq_record.id}_{row['database']}_{row['query_file_name']}_full_contig"
            description = "full contig extracted"
            return SeqIO.SeqRecord(seq_record.seq, id=header_id, description=description)

    if not found:
        msg = f"Sequence ID {row['sseqid']} not found in {original_fasta}"
        print(f"\n\033[91m{msg}\033[0m")
        logger.warning(msg)
    return None

def extract_contigs_from_csv(csv_path, fasta_input, output_fasta, evalue=1e-5,
                             min_perc=90.0, min_cov=75.0, *, keep_temp_files=False, logger=None,
                             threads):
    df = pd.read_csv(csv_path)
    contigs = []
    extracted_contigs = set()
    lock = Lock()

    #Collect all FASTA files (file_handler)
    try:
        fasta_files, temp_dir = collect_fasta_files(fasta_input)
    except ValueError as e:
        print(f"\n\033[91mError: {e}\033[0m")
        logger.error(str(e))
        return

    #Build a lookup map for quick database to FASTA resolution
    fasta_map = {}
    for f in fasta_files:
        base = os.path.splitext(os.path.basename(f))[0]
        fasta_map[base] = f

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(
                process_contig_entry,
                row, fasta_map, evalue, min_perc, min_cov,
                extracted_contigs, lock, logger=logger
            )
            for _, row in df.iterrows()
        ]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                contigs.append(result)

    SeqIO.write(contigs, output_fasta, "fasta")
    print(f"\033[92mExtracted {len(contigs)} unique contigs to {output_fasta}\033[0m")
    logger.info(f"Extracted {len(contigs)} unique contigs to {output_fasta}")

    #Clean up any extracted archive temp files (/tmp/)
    cleanup_temp_dir(temp_dir, keep=keep_temp_files, logger=logger)

def run_contigs(args):
    logger = logging.getLogger("extract-contig")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if not any(isinstance(h, logging.FileHandler) and h.baseFilename.endswith("seqforge_extract-contig.log")
               for h in logger.handlers):
        file_handler = logging.FileHandler("seqforge_extract-contig.log", mode='a')
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    start_time = datetime.now()
    print(f"Contig extraction started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Contig extraction started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    extract_contigs_from_csv(
        csv_path=args.csv_path,
        fasta_input=args.fasta_directory,
        output_fasta=args.output_fasta,
        evalue=args.evalue,
        min_perc=args.min_perc,
        min_cov=args.min_cov,
        keep_temp_files=args.keep_temp_files,
        logger=logger,
        threads=args.threads if args.threads else 4
    )

    end_time = datetime.now()
    print(f"Contig extraction completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Contig extraction completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    msg = f"Total time: {str(end_time - start_time)}"
    logger.info(msg)
    print(msg)
    

# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import os
from datetime import datetime
import logging
import pandas as pd
from Bio import SeqIO
import concurrent.futures
from threading import Lock

def find_fasta_file(basename, fasta_input):
    extensions = {'fa', 'fas', 'fasta', 'fna', 'faa'}
    
    if os.path.isfile(fasta_input):
        name_part, extension = os.path.splitext(os.path.basename(fasta_input))
        extension = extension.lstrip('.').lower()
        if extension in extensions and (
            name_part == basename or
            os.path.basename(fasta_input).startswith(basename + '_') or
            os.path.basename(fasta_input).startswith(basename + '.')
        ):
            return fasta_input
        return None

    for file in os.listdir(fasta_input):
        name_part, extension = os.path.splitext(file)
        extension = extension.lstrip('.').lower()
        if extension in extensions:
            if name_part == basename or file.startswith(basename + '_') or file.startswith(basename + '.'):
                return os.path.join(fasta_input, file)

    return None

def process_contig_entry(row, fasta_input, evalue, min_perc, min_cov, extracted_contigs, lock, *, logger):
    if float(row['evalue']) > evalue or float(row['pident']) < min_perc or float(row['query_coverage']) < min_cov:
        msg = f"Skipping sequence ({row['database']} query: {row['query_file_name']} pident: {row['pident']} query coverage: {row['query_coverage']} evalue: {row['evalue']}) due to filtering thresholds"
        print(f"\n\033[91m{msg}\033[0m")
        logger.info(msg)
        return None
    
    contig_key = (row['database'], row['sseqid'])

    with lock:
        if contig_key in extracted_contigs:
            return None
        extracted_contigs.add(contig_key)
    
    original_fasta = find_fasta_file(row['database'], fasta_input)
    if original_fasta is None:
        msg = f"No matching FASTA file found for {row['database']}"
        print(f"\033[91m{msg}\033[0m")
        logger.warning(msg)
        return None
    
    found = False
    for seq_record in SeqIO.parse(original_fasta, "fasta"):
        seq_id = seq_record.id.split()[0]
        if seq_record.id == row['sseqid']:
            found = True
            extracted_contigs.add(contig_key) # Mark this contig as extracted
            header_id = f"{seq_record.id}_{row['database']}_{row['query_file_name']}_full_contig"
            description = "full contig extracted"
            return SeqIO.SeqRecord(seq_record.seq, id=header_id, description=description)
        
    if not found:
        msg = f"Sequence ID {row['sseqid']} not found in {original_fasta}Sequence ID {row['sseqid']} not found in {original_fasta}"
        print(f"\n\033[91m{msg}\033[0m")
        logger.warning(msg)
    return None
    
def extract_contigs_from_csv(csv_path, fasta_input, output_fasta, evalue=1e-5, min_perc=90.0, min_cov=75.0,
                             *, logger):
    df = pd.read_csv(csv_path)
    contigs = []
    extracted_contigs = set() # Set to keep track of contigs that have already been extracted
    lock = Lock()

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(
                process_contig_entry, 
                row, fasta_input, evalue, min_perc, 
                min_cov, extracted_contigs, lock, 
                logger=logger
            )
            for index, row in df.iterrows()
        ]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                contigs.append(result)

    SeqIO.write(contigs, output_fasta, "fasta")
    print(f"\n\033[92mExtracted {len(contigs)} unique contigs to {output_fasta}\033[0m")
    logger.info(f"Extracted {len(contigs)} unique contigs to {output_fasta}")

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
        logger=logger
    )

    end_time = datetime.now()
    print(f"Contig extraction completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Contig extraction completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    msg = f"Total time: {str(end_time - start_time)}"
    logger.info(msg)
    print(msg)
    
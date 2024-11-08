import os
import time
import pandas as pd
from Bio import SeqIO
import concurrent.futures
from threading import Lock

"""
This program is designed to extract entire contigs based on a multiBLAST Query.

Author: Elijah R. Bring Horvath, PhD
"""

def find_fasta_file(basename, fasta_dir):
    extensions = {'fa', 'fna', 'fas', 'fasta', 'faa'}
    for file in os.listdir(fasta_dir):
        name_part, extension = os.path.splitext(file)
        extension = extension.lstrip('.').lower()
        if extension in extensions:
            if name_part == basename or file.startswith(basename + '_') or file.startswith(basename + '.'):
                return os.path.join(fasta_dir, file)
    return None

def process_contig_entry(row, fasta_dir, min_evalue, min_perc, min_cov, extracted_contigs, lock):
    if float(row['evalue']) > min_evalue or float(row['pident']) < min_perc or float(row['query_coverage']) < min_cov:
        print(f"\n \033[91mSkipping sequence ({row['database']} query: {row['query_file_name']} pident: {row['pident']} query coverage: {row['query_coverage']} evalue: {row['evalue']}) due to filtering thresholds\033[0m")
        return None
    
    contig_key = (row['database'], row['sseqid'])

    with lock:
        if contig_key in extracted_contigs:
            return None
        extracted_contigs.add(contig_key)
    
    original_fasta = find_fasta_file(row['database'], fasta_dir)
    if original_fasta is None:
        print(f"\n \033[91mNo matching FASTA file found for {row['database']}\033[0m")
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
        print(f"\n \033[91mSequence ID {row['sseqid']} not found in {original_fasta}\033[0m")
    return None
    
def extract_contigs_from_csv(csv_path, fasta_dir, output_fasta, min_evalue=1e-5, min_perc=90.0, min_cov=75.0):
    df = pd.read_csv(csv_path)
    contigs = []
    extracted_contigs = set() # Set to keep track of contigs that have already been extracted
    lock = Lock()

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(process_contig_entry, row, fasta_dir, min_evalue, min_perc, min_cov, extracted_contigs, lock)
            for index, row in df.iterrows()
        ]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                contigs.append(result)

    SeqIO.write(contigs, output_fasta, "fasta")
    print(f"\n \033[92mExtracted {len(contigs)} unique contigs to {output_fasta}\033[0m")

def run_contigs(args):
    start_time = time.time()

    extract_contigs_from_csv(
        csv_path=args.csv_path,
        fasta_dir=args.fasta_directory,
        output_fasta=args.output_fasta,
        min_evalue=args.min_evalue,
        min_perc=args.min_perc,
        min_cov=args.min_cov
    )

    end_time = time.time()
    print(f" Total runtime: {end_time-start_time:.2f} seconds")
    

import os
import time
import pandas as pd
from Bio import SeqIO
import concurrent.futures

"""
This program is designed to extract sequences based on a multiBLAST query. Translation of sequences is optional, however care
    should be used when translating extracted sequences, as BLAST results may not always contain a full CDS. To allow for this, when 
    the --translate argument is called, extracted sequences will be trimmed to only include complete codons, which may affect interpretation
    of results. 
"""
#Process input FASTA files
def find_fasta_file(basename, fasta_dir):
    
    extensions = {'fa', 'fas', 'fasta', 'fna', 'faa'}
    for file in os.listdir(fasta_dir):
        name_part, extension = os.path.splitext(file)
        extension = extension.lstrip('.').lower()
        if extension in extensions:
            if name_part == basename or file.startswith(basename + '_') or file.startswith(basename + '.'):
                return os.path.join(fasta_dir, file)
    return None

def process_sequence_entry(row, fasta_dir, translate,
                           min_evalue, min_perc, min_cov,
                           up, down):
    """Extract aligned region ± up/down bp, apply filters & optional translation."""
    # filter by BLAST thresholds
    if (float(row['evalue']) > min_evalue or
        float(row['pident']) < min_perc or
        float(row['query_coverage']) < min_cov):
        print(f"\n \033[91mSkipping sequence ({row['database']} query: {row['query_file_name']} "
              f"pident: {row['pident']} cov: {row['query_coverage']} evalue: {row['evalue']})\033[0m")
        return None

    # find the correct FASTA file
    original_fasta = find_fasta_file(row['database'], fasta_dir)
    if original_fasta is None:
        print(f"\n \033[91mNo matching FASTA file found for {row['database']}\033[0m")
        return None

    # parse and extract region
    for seq_record in SeqIO.parse(original_fasta, "fasta"):
        if seq_record.id != row['sseqid']:
            continue

        sstart, send = int(row['sstart']), int(row['send'])
        low, high = min(sstart, send), max(sstart, send)
        contig_len = len(seq_record.seq)

        ideal_start = low - up
        ideal_end = high + down

        # extend but clamp to contig boundaries
        region_start = max(1, low  - up)
        region_end   = min(contig_len, high + down)

        if ideal_start < 1:
            print(f" \033[93mWarning: Requested upstream extension of {up} bp for "
                  f"{seq_record.id} truncated at contig start (only {low-1} bp available\033[0m)")
        if ideal_end > contig_len:
            print(f" \033[93mWarning: Requested downstream extension of {down} bp for "
                  f"{seq_record.id} truncated at contig end (only {contig_len-high} bp available)\033[0m")

        subseq = seq_record.seq[region_start-1:region_end]
        # reverse-complement if hit on negative strand
        if sstart > send:
            subseq = subseq.reverse_complement()

        # optional translation (trim to codon)
        if translate:
            trim_len = len(subseq) - (len(subseq) % 3)
            subseq = subseq[:trim_len].translate()

        header_id   = f"{seq_record.id}_{row['database']}_{row['query_file_name']}_region"
        description = "translated" if translate else "nucleotide"
        return SeqIO.SeqRecord(subseq, id=header_id, description=description)

    # if we get here, no matching record ID was found
    print(f"\n \033[91mSequence ID {row['sseqid']} not found in {original_fasta}\033[0m")
    return None

def extract_sequences_from_csv(csv_path, fasta_dir, output_fasta,
                               translate=False, min_evalue=1e-5,
                               min_perc=90.0, min_cov=75.0,
                               up=0, down=0):
    """Read BLAST CSV, extract regions in parallel, write to output FASTA."""
    df = pd.read_csv(csv_path)
    sequences = []

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(
                process_sequence_entry,
                row, fasta_dir, translate,
                min_evalue, min_perc, min_cov,
                up, down
            )
            for _, row in df.iterrows()
        ]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                sequences.append(result)

    SeqIO.write(sequences, output_fasta, "fasta")
    print(f"\n \033[92mExtracted {len(sequences)} sequences to {output_fasta}\033[0m\n")

def run(args):
    """Entry point for 'extract' sub-command."""
    start_time = time.time()
    extract_sequences_from_csv(
        csv_path     = args.csv_path,
        fasta_dir    = args.fasta_directory,
        output_fasta = args.output_fasta,
        translate    = args.translate,
        min_evalue   = args.min_evalue,
        min_perc     = args.min_perc,
        min_cov      = args.min_cov,
        up           = args.up,
        down         = args.down,
    )
    end_time = time.time()
    print(f" Total runtime: {end_time - start_time:.2f} seconds")
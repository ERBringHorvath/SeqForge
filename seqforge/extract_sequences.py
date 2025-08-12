# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import os
import logging
from datetime import datetime
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import concurrent.futures

from utils.file_handler import collect_fasta_files, cleanup_temp_dir

NUC_CHARS = set("ACGTRYSWKMBDHVNU-.") #IUPAC nucleotides + gap/dot

def is_protein_fasta(path, max_records=25, max_chars=20000):
    seen = set()
    total = 0
    try:
        for i, rec in enumerate(SeqIO.parse(path, "fasta"), start=1):
            s = str(rec.seq).upper()
            for ch in s:
                if ch.isalpha() or ch in "-.*":
                    seen.add(ch)
                    total += 1
                    if ch not in NUC_CHARS:
                        return True #definitely protein
                    if total >= max_chars:
                        return False
            if i >= max_records:
                break
    except Exception:
        return False
    #If we only see NT chars, call it NT
    return False

def process_sequence_entry(row, fasta_map, translate,
                           evalue, min_perc, min_cov,
                           up, down, *, logger):

    if (float(row['evalue']) > evalue or
        float(row['pident']) < min_perc or
        float(row['query_coverage']) < min_cov):
        print(f"\033[91mSkipping sequence ({row['database']} query: {row['query_file_name']} "
              f"pident: {row['pident']} cov: {row['query_coverage']} evalue: {row['evalue']})\033[0m")
        return None

    original_fasta = fasta_map.get(row['database'])
    if not original_fasta:
        msg = f"No matching FASTA file found for {row['database']}"
        print(f"\033[91m{msg}\033[0m")
        logger.warning(msg)
        return None

    #parse and extract region
    with open(original_fasta, "r") as handle:
        for seq_record in SeqIO.parse(handle, "fasta"):
            if seq_record.id != row['sseqid']:
                continue

            sstart, send = int(row['sstart']), int(row['send'])
            low, high = min(sstart, send), max(sstart, send)
            contig_len = len(seq_record.seq)

            ideal_start = low - up
            ideal_end = high + down

            #clamp to contig boundaries
            region_start = max(1, low  - up)
            region_end   = min(contig_len, high + down)

            if ideal_start < 1:
                print(f"\033[93mWarning: Requested upstream extension of {up} bp for "
                    f"{seq_record.id} truncated at contig start (only {low-1} bp available)\033[0m")
                logger.warning(f"Warning: Requested upstream extension of {up} bp for "
                            f"{seq_record.id} truncated at contig start (only {low-1} bp available)")
            if ideal_end > contig_len:
                print(f"\033[93mWarning: Requested downstream extension of {down} bp for "
                    f"{seq_record.id} truncated at contig end (only {contig_len-high} bp available)\033[0m")
                logger.warning(f"Warning: Requested downstream extension of {down} bp for "
                            f"{seq_record.id} truncated at contig end (only {contig_len-high} bp available)")

            subseq = seq_record.seq[region_start-1:region_end]
            #reverse-complement if hit on negative strand
            if sstart > send:
                subseq = subseq.reverse_complement()

            if translate:
                trim_len = len(subseq) - (len(subseq) % 3)
                subseq = subseq[:trim_len].translate()

            header_id   = f"{seq_record.id}_{row['database']}_{row['query_file_name']}_region"
            description = "translated" if translate else "nucleotide"
            return SeqRecord(subseq, id=header_id, description=description)

    msg = f"Sequence ID {row['sseqid']} not found in {original_fasta}"
    print(f"\033[91m{msg}\033[0m")
    logger.warning(msg)
    return None

def extract_sequences_from_csv(csv_path, fasta_input, output_fasta,
                               translate=False, evalue=1e-5,
                               min_perc=90.0, min_cov=75.0,
                               up=0, down=0, *, keep_temp_files=False, logger=None,
                               threads=4, temp_dir_base=None):
    
    df = pd.read_csv(csv_path)
    sequences = []

    #Collect all FASTA files (file_handler)
    try:
        fasta_files, temp_dir = collect_fasta_files(fasta_input, temp_dir_base=temp_dir_base)
    except ValueError as e:
        print(f"\n\033[91mError: {e}\033[0m")
        logger.error(str(e))
        return
    
    if translate:
        offending = []
        for f in fasta_files:
            name = os.path.basename(f).lower()
            if name.endswith((".faa", ".faa.gz")) or is_protein_fasta(f):
                offending.append(f)
        if offending:
            msg = {
                "Amino acid FASTA(s) detected; --translate expects nucleotide sequences\n"
                "Offending files:\n " + "\n ".join(offending[:10]) + ("..." if len(offending) > 10 else "")
            }
            print(f"\033[91m{msg}\033[0m")
            if logger:
                logger.warning(msg)
            cleanup_temp_dir(temp_dir, keep=keep_temp_files, logger=logger)
            return

    fasta_map = {}
    for f in fasta_files:
        base = os.path.splitext(os.path.basename(f))[0]
        fasta_map[base] = f

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        # up = int(up)
        # down = int(down)
        futures = [
            executor.submit(
                process_sequence_entry,
                row, fasta_map, translate,
                evalue, min_perc, min_cov,
                int(up), int(down), logger=logger
            )
            for _, row in df.iterrows()
        ]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                sequences.append(result)

    SeqIO.write(sequences, output_fasta, "fasta")
    print(f"\033[92mExtracted {len(sequences)} sequences to {output_fasta}\033[0m\n")
    logger.info(f"Extracted {len(sequences)} sequences to {output_fasta}")

    #Clean up temp directory from archive extraction if archive input
    cleanup_temp_dir(temp_dir, keep=keep_temp_files, logger=logger)

def run(args):
    logger = logging.getLogger("extract")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if not any(isinstance(h, logging.FileHandler) and h.baseFilename.endswith("seqforge_extract.log")
               for h in logger.handlers):
        file_handler = logging.FileHandler("seqforge_extract.log", mode='a')
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    start_time = datetime.now()
    print(f"Extraction started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    logger.info(f"Extraction started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    if args.translate:
        if args.up:
            msg = "Error: Basepair padding does not currently support translation"
            print(f"\n\033[91m{msg}\033[0m")
            logger.warning(msg)
            return
        if args.down:
            msg = "Error: Basepair padding does not currently support translation"
            print(f"\n\033[91m{msg}\033[0m")
            logger.warning(msg)
            return

    extract_sequences_from_csv(
        csv_path=args.csv_path,
        fasta_input=args.fasta_directory,
        output_fasta=args.output_fasta,
        translate=args.translate,
        evalue=args.evalue,
        min_perc=args.min_perc,
        min_cov=args.min_cov,
        up=int(args.up),
        down=int(args.down),
        keep_temp_files=getattr(args, 'keep_temp_files', False),
        temp_dir_base=getattr(args, 'temp_dir', None),
        logger=logger,
        threads=args.threads if args.threads else 4
    )

    end_time = datetime.now()
    print(f"Extraction completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Extraction completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    msg = f"Total time: {str(end_time - start_time)}"
    logger.info(msg)
    print(msg)

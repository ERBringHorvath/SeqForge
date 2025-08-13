# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import os 
import logging 
import gzip
from io import TextIOWrapper
from Bio import SeqIO
from math import ceil
from datetime import datetime

from .utils.progress import ProgressHandler

def run_split(args):
    logger = logging.getLogger("split-fasta")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if not any(isinstance(h, logging.FileHandler) and h.baseFilename.endswith("seqforge_split-fasta.log")
            for h in logger.handlers):
        file_handler = logging.FileHandler("seqforge_split-fasta.log", mode='a')
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    start_time = datetime.now()
    print(f"Split-FASTA started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Split-FASTA started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    #arg wrangling
    input_file = args.input
    output_dir = args.output_dir
    fragment_size = args.fragment
    compress = args.compress

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    allowed_extensions = ['.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa']
    input_basename = os.path.basename(input_file)
    base_name, _ = os.path.splitext(input_basename)

    _, file_extension = os.path.splitext(input_file)
    if file_extension.lower() not in allowed_extensions:
        logging.warning(f"\n \033[91mSkipped {input_file} due to incorrect file extension.\033[0m")
        return

    try:
        records = list(SeqIO.parse(input_file, "fasta"))
    except Exception as e:
        logging.error(f"\n \033[91mError parsing {input_file}: {e}\033[0m")
        return

    total_records = len(records)

    if fragment_size:
        num_chunks = ceil(total_records / fragment_size) #round to closest int
        digits = len(str(num_chunks))
        msg = f"Splitting {total_records} sequences into {num_chunks} file(s) of ~{fragment_size} each"
        print(f"\n\033[94m{msg}\033[0m")
        logger.info(msg)

        progress = ProgressHandler(total=total_records, prefix="Fragmenting", mode='bar')

        for i in range(num_chunks):
            chunk_records = records[i * fragment_size : (i + 1) * fragment_size]
            chunk_num = str(i + 1).zfill(digits)
            ext = ".fa.gz" if compress else ".fa"
            out_name = f"{base_name}_Frag{fragment_size}_{chunk_num}{ext}"
            out_path = os.path.join(output_dir, out_name)
            try:
                if compress:
                    with gzip.open(out_path, "wt") as f:
                        SeqIO.write(chunk_records, f, "fasta")
                else:
                    SeqIO.write(chunk_records, out_path, "fasta")
                progress.update(len(chunk_records))
            except Exception as e:
                logging.error(f"\n \033[91mError writing chunk {chunk_num}: {e}\033[0m")
        
        progress.finish()
        print()
    else:
        msg = "Splitting into individual FASTA files..."
        print(f"\n\033[94m{msg}\033[0m")
        logger.info(msg)

        progress = ProgressHandler(total=total_records, prefix="Processing", mode='bar')

        for record in records:
            out_name = f"{record.id}.fasta.gz" if compress else f"{record.id}.fasta"
            out_path = os.path.join(output_dir, out_name)
            try:
                if compress:
                    with gzip.open(out_path, "wt") as f:
                        SeqIO.write([record], f, "fasta")
                else:
                    SeqIO.write(record, out_path, "fasta")
                progress.update(1, current_item=record.id)
            except Exception as e:
                logging.error(f"\n \033[91mError writing {record.id}: {e}\033[0m")

        progress.finish()
        print()

    if args.compress:
        print("\n\033[93mCompressed output files (.fasta.gz)\033[0m\n")
        logger.info("Compressed output files (.fasta.gz)")

    end_time = datetime.now()
    print(f"Split-FASTA completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Split-FASTA completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    msg = f"Total runtime: {str(end_time - start_time)}"
    logger.info(msg)
    print(msg)

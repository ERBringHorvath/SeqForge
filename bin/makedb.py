# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import os
import subprocess
import logging
import gzip
from concurrent.futures import ProcessPoolExecutor

logger = logging.getLogger("makedb")
logger.setLevel(logging.INFO)

if not logger.handlers:
    file_handler = logging.FileHandler("seqforge_makedb.log", mode='w')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

def makeblastdb_single(input_file, output_dir):
    filename = os.path.basename(input_file)
    basename = os.path.splitext(filename)[0]

    if filename.endswith(".gz"):
        basename = os.path.splitext(basename)[0]  # Remove .gz and .fasta

    output_file = os.path.join(output_dir, basename)

    if input_file.endswith((".fasta", ".fna", ".ffn", ".fa", ".fas", ".fasta.gz", ".fna.gz", ".ffn.gz", ".fa.gz", ".fas.gz")):
        dbtype = "nucl"
    elif input_file.endswith((".faa", ".faa.gz")):
        dbtype = "prot"
    else:
        print(f"\n\033[91mWarning: Unrecognized extension for {filename}. Skipping\033[0m")
        logging.warning(f"Skipped unrecognized file: {input_file}")
        return

    if input_file.endswith(".gz"):
        cmd = (
            f"gunzip -c \"{input_file}\" | "
            f"makeblastdb -in - -out \"{output_file}\" -dbtype {dbtype} -title \"{basename}\""
    )
    else:
        cmd = f"makeblastdb -in \"{input_file}\" -out \"{output_file}\" -dbtype {dbtype}"

    # Suppress stdout, capture stderr
    result = subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        logger.error(f"Failed to create DB for {input_file}: {result.stderr.strip()}")
    else:
        logger.info(f"Successfully created DB for {input_file}")

def make_blast_db(args):
    input_path = args.fasta_directory
    output_dir = args.out
    threads = args.threads if hasattr(args, 'threads') and args.threads else 4

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    valid_exts = (".fasta", ".faa", ".fna", ".ffn", ".fa", ".fas",
                  ".fasta.gz", ".faa.gz", ".fna.gz", ".ffn.gz", ".fa.gz", ".fas.gz")

    files_to_process = []
    if os.path.isdir(input_path):
        files_to_process = [
            os.path.join(input_path, f)
            for f in os.listdir(input_path)
            if f.endswith(valid_exts)
        ]
    elif os.path.isfile(input_path) and input_path.endswith(valid_exts):
        files_to_process = [input_path]
    else:
        print(f"\n\033[91mError: No valid FASTA file(s) found at {input_path}\033[0m")
        logger.error(f"No valid FASTA files found at input path: {input_path}")
        return

    with ProcessPoolExecutor(max_workers=threads) as executor:
        executor.map(makeblastdb_single, files_to_process, [output_dir]*len(files_to_process))

    print("\n\033[92mDatabases successfully created\033[0m")
    logger.info("Finished creating BLAST databases.")

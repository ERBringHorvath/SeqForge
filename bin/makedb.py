# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import os
import subprocess
import logging
from datetime import datetime
import re
from concurrent.futures import ProcessPoolExecutor

# Valid FASTA extensions
FASTA_EXTENSIONS = (".fasta", ".faa", ".fna", ".ffn", ".fa", ".fas",
                    ".fasta.gz", ".faa.gz", ".fna.gz", ".ffn.gz", ".fa.gz", ".fas.gz")

def sanitize_filename(filename):
    name, ext = os.path.splitext(filename)
    name = name.replace('.', '_')
    name = re.sub(r'[-:;\s]', '_', name)
    name = re.sub(r"[()'`\"‘’“”]", '', name)
    name = re.sub(r'_+', '_', name)
    return name + ext

def contains_special_characters(filenames):
    return any(sanitize_filename(f) != f for f in filenames)

def sanitize_files_in_place(filepaths):
    for path in filepaths:
        dir_path = os.path.dirname(path)
        old_name = os.path.basename(path)
        new_name = sanitize_filename(old_name)
        if old_name != new_name:
            os.rename(os.path.join(dir_path, old_name), os.path.join(dir_path, new_name))

def makeblastdb_single(input_file, output_dir):
    filename = os.path.basename(input_file)
    basename = os.path.splitext(filename)[0]

    if filename.endswith(".gz"):
        basename = os.path.splitext(basename)[0]

    output_file = os.path.join(output_dir, basename)

    if input_file.endswith((".fasta", ".fna", ".ffn", ".fa", ".fas", ".fasta.gz", ".fna.gz", ".ffn.gz", ".fa.gz", ".fas.gz")):
        dbtype = "nucl"
    elif input_file.endswith((".faa", ".faa.gz")):
        dbtype = "prot"
    else:
        print(f"\n\033[91mWarning: Unrecognized extension for {filename}. Skipping\033[0m")
        logger.warning(f"Skipped unrecognized file: {input_file}")
        return

    if input_file.endswith(".gz"):
        cmd = f"gunzip -c \"{input_file}\" | makeblastdb -in - -out \"{output_file}\" -dbtype {dbtype} -title \"{basename}\""
    else:
        cmd = f"makeblastdb -in \"{input_file}\" -out \"{output_file}\" -dbtype {dbtype}"

    result = subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        logger.error(f"Failed to create DB for {input_file}: {result.stderr.strip()}")
    else:
        logger.info(f"Successfully created DB for {input_file}")

def run_make_blast_db(args):

    logger = logging.getLogger("makedb")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if not any(isinstance(h, logging.FileHandler) and h.baseFilename.endswith("seqforge_makedb.log")
            for h in logger.handlers):
        file_handler = logging.FileHandler("seqforge_makedb.log", mode='a')
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    logger.info("Starting Database Creation")
    start_time = datetime.now()

    input_path = args.fasta_directory
    output_dir = args.output_dir
    threads = args.threads if args.threads else 4
    sanitize_flag = getattr(args, 'sanitize', False)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Gather files
    files_to_process = []
    if os.path.isdir(input_path):
        files_to_process = [
            os.path.join(input_path, f)
            for f in os.listdir(input_path)
            if f.endswith(FASTA_EXTENSIONS)
        ]
    elif os.path.isfile(input_path) and input_path.endswith(FASTA_EXTENSIONS):
        files_to_process = [input_path]
    else:
        print(f"\n\033[91mError: No valid FASTA file(s) found at {input_path}\033[0m")
        logger.error("No valid FASTA files found")
        return

    filenames_only = [os.path.basename(f) for f in files_to_process]

    if contains_special_characters(filenames_only):
        if sanitize_flag:
            sanitize_files_in_place(files_to_process)
            # Re-scan files
            files_to_process = [
                os.path.join(input_path, sanitize_filename(os.path.basename(f)))
                for f in files_to_process
            ]
            print("\n\033[93mSanitized filenames in place\033[0m")
        else:
            print("\n\033[91mError: Filenames contain special characters\033[0m")
            print("\033[91mPlease rerun with --sanitize (permanently removes special characters from file names)" 
                  " or first run the 'seqforge sanitize' module\033[0m")
            print("\033[93mRun 'seqforge sanitize --help' for more information\033[0m")
            logger.error("Filenames contain special characters. Aborting")
            return

    with ProcessPoolExecutor(max_workers=threads) as executor:
        executor.map(makeblastdb_single, files_to_process, [output_dir]*len(files_to_process))

    print("\n\033[92mDatabases successfully created\033[0m")
    logger.info("Finished creating BLAST databases")
    end_time = datetime.now()
    print(f"Completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total runtime: {str(end_time - start_time)}")

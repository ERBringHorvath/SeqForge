# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import os
import shutil
import csv
import gzip
import re
import random
import logging
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO
from Bio.Seq import Seq

def generate_scramble(length):
    return ''.join(random.choices('ACGT', k=length))

def load_sequences(sequence_file, output_dir, logger):
    valid_sequences = set()
    skipped_sequences = []
    duplicate_sequences = []
    total_lines = 0

    open_func = gzip.open if sequence_file.endswith((".gz", ".gzip")) else open

    def add_sequence(seq):
        seq = seq.upper().strip()
        if not seq or not re.fullmatch(r"[ACGT]+", seq):
            skipped_sequences.append(seq)
            return
        rc = str(Seq(seq).reverse_complement())
        for s in (seq, rc):
            if s in valid_sequences:
                duplicate_sequences.append(s)
            valid_sequences.add(s)

    try:
        # Try to parse as FASTA
        with open_func(sequence_file, 'rt') as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            if records:
                print("\033[94mFASTA format detected\033[0m")
                logger.info("FASTA format detected")
                for record in records:
                    total_lines += 1
                    add_sequence(str(record.seq))
                return sorted(valid_sequences)
    except Exception:
        pass  # Fallback to text mode below

    # If FASTA parse failed, treat as plain text
    print("\n\033[93mPlain text format detected\033[0m")
    logger.info("Plain text format detected")
    with open_func(sequence_file, 'rt') as f:
        for line in f:
            total_lines += 1
            add_sequence(line.strip())

    if skipped_sequences:
        skipped_path = os.path.join(output_dir, "skipped_sequences.txt")
        with open(skipped_path, "w") as f:
            for entry in skipped_sequences:
                f.write(entry + "\n")
        msg = f"Skipped {len(skipped_sequences)} invalid sequences. Logged to: {skipped_path}"
        print(f"\n\033[91m{msg}\033[0m")
        logger.info(msg)

    if duplicate_sequences:
        msg = f"{len(duplicate_sequences)} duplicate (or reverse complement) sequences were ignored"
        print(f"\n\033[91m{msg}\033[0m")
        logger.info(msg)
        print("\033[91mDuplicates or Reverse Compliments:\033[0m")
        logger.info("Duplicates or Reverse Compliments:")
        for dup in sorted(set(duplicate_sequences)):
            print(f" - {dup}")
            logger.info(f" - {dup}")

    return sorted(valid_sequences)

def _mask_single_file(args):
    filename, input_dir, output_dir, sequences, verbose, dash, scramble, dry_run = args

    input_path = os.path.join(input_dir, filename)
    output_path = os.path.join(output_dir, filename)
    replacement_char = '-' if dash else 'N'

    modified = False
    matched_sequences = set()
    records = []

    try:
        for record in SeqIO.parse(input_path, "fasta"):
            sequence = str(record.seq).upper()

            for seq in sequences:
                if seq in sequence:
                    if scramble:
                        if len(seq) < 3:
                            masked_seq = generate_scramble(len(seq))
                        else:
                            core = generate_scramble(len(seq) - 2)
                            masked_seq = f"-{core}-"
                    else:
                        masked_seq = replacement_char * len(seq)

                    sequence = sequence.replace(seq, masked_seq)
                    matched_sequences.add(seq)
                    modified = True

            record.seq = Seq(sequence)
            records.append(record)

        if modified:
            if dry_run:
                print(f"[DRY RUN] Would mask: {filename}")
            else:
                SeqIO.write(records, output_path, "fasta")
        else:
            if not dry_run:
                shutil.copy(input_path, output_path)

        return {
            "filename": filename,
            "modified": modified,
            "matched_sequences": sorted(matched_sequences) if verbose else None
        }

    except Exception as e:
        return {"filename": filename, "error": str(e)}

def mask_sequences(input_dir, output_dir, sequences, logger, verbose=False, dash=False, scramble=False, threads=1, dry_run=False):
    os.makedirs(output_dir, exist_ok=True)
    fasta_extensions = {".fasta", ".fa", ".fna", ".fas", ".ffn"}
    fasta_files = [f for f in os.listdir(input_dir) if any(f.lower().endswith(ext) for ext in fasta_extensions)]

    summary_log = []
    masked_files = []

    print(f"\n\033[92mProcessing {len(fasta_files)} FASTA files using {threads} thread(s)...\033[0m")
    logger.info(f"Processing {len(fasta_files)} FASTA files using {threads} thread(s)...")

    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(_mask_single_file, (f, input_dir, output_dir, sequences, verbose, dash, scramble, dry_run))
            for f in fasta_files
        ]

        for future in as_completed(futures):
            result = future.result()
            if result.get("error"):
                print(f"\033[91mError in {result['filename']}: {result['error']}\033[0m")
                logger.warning(f"Error in {result['filename']}: {result['error']}")
                continue
            if result["modified"]:
                masked_files.append(result["filename"])
            if verbose and result["modified"]:
                summary_log.append({
                    "filename": result["filename"],
                    "matched_sequences": ",".join(result["matched_sequences"])
                })
    if not dry_run:
        print(f"\n\033[93m{'[DRY RUN] ' if dry_run else ''}Summary:\033[0m")
        logger.info(f"{'[DRY RUN] ' if dry_run else ''}")
        print(f"\033[92mTotal FASTA files checked: {len(fasta_files)}\033[0m")
        print(f"\033[92mFiles with at least one sequence match: {len(masked_files)}\033[0m")
        if masked_files:
            print("\033[92mFiles with matches:\033[0m")
            logger.info("Files with matches:")
            for f in sorted(masked_files):
                print(f" - {f}")
                logger.info(f" - {f}")
        else:
            print("\n\033[91mNo files contained sequence matches\033[0m")

    if verbose and summary_log:
        csv_path = os.path.join(output_dir, "masked_summary.csv")
        with open(csv_path, "w", newline="") as csvfile:
            fieldnames = ["filename", "matched_sequences"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for row in summary_log:
                writer.writerow(row)
            print(f"\n\033[95mVerbose log saved to: {csv_path}\n\033[0m")
            logger.info(f"Verbose log saved to : {csv_path}")

def run_mask(args):
    logger = logging.getLogger("mask")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if not any(isinstance(h, logging.FileHandler) and h.baseFilename.endswith("seqforge_mask.log")
            for h in logger.handlers):
        file_handler = logging.FileHandler("seqforge_mask.log", mode='a')
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    start_time = datetime.now()
    print(f"Masking started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Masking started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    os.makedirs(args.output_dir, exist_ok=True)
    sequences = load_sequences(args.sequence_file, args.output_dir, logger)
    mask_sequences(
        args.input_dir,
        args.output_dir,
        sequences,
        logger,
        verbose=args.verbose,
        dash=args.dash,
        scramble=args.scramble,
        threads=args.threads,
        dry_run=args.dry_run
    )

    end_time = datetime.now()
    print(f"Masking completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Masking completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    msg = f"Total runtime: {str(end_time - start_time)}"
    logger.info(msg)
    print(msg)

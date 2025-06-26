# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import os
import shutil
import csv
import gzip
import re
import random
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO
from Bio.Seq import Seq

def generate_scramble(length):
    return ''.join(random.choices('ACGT', k=length))

def load_sequences(sequence_file, output_dir):
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
                for record in records:
                    total_lines += 1
                    add_sequence(str(record.seq))
                return sorted(valid_sequences)
    except Exception:
        pass  # Fallback to text mode below

    # If FASTA parse failed, treat as plain text
    print("\n\033[93mPlain text format detected\033[0m")
    with open_func(sequence_file, 'rt') as f:
        for line in f:
            total_lines += 1
            add_sequence(line.strip())

    if skipped_sequences:
        skipped_path = os.path.join(output_dir, "skipped_sequences.txt")
        with open(skipped_path, "w") as f:
            for entry in skipped_sequences:
                f.write(entry + "\n")
        print(f"\n033\91mSkipped {len(skipped_sequences)} invalid sequences. Logged to: {skipped_path}\033[0m")

    if duplicate_sequences:
        print(f"\n\033[91m{len(duplicate_sequences)} duplicate (or reverse complement) sequences were ignored\033[0m")
        print("\033[91mDuplicates or Reverse Compliments:\033[0m")
        for dup in sorted(set(duplicate_sequences)):
            print(f" - {dup}")

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

def mask_sequences(input_dir, output_dir, sequences, verbose=False, dash=False, scramble=False, threads=1, dry_run=False):
    os.makedirs(output_dir, exist_ok=True)
    fasta_extensions = {".fasta", ".fa", ".fna", ".fas", ".ffn"}
    fasta_files = [f for f in os.listdir(input_dir) if any(f.lower().endswith(ext) for ext in fasta_extensions)]

    summary_log = []
    masked_files = []

    print(f"\n\033[92mProcessing {len(fasta_files)} FASTA files using {threads} thread(s)...\033[0m")

    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(_mask_single_file, (f, input_dir, output_dir, sequences, verbose, dash, scramble, dry_run))
            for f in fasta_files
        ]

        for future in as_completed(futures):
            result = future.result()
            if result.get("error"):
                print(f"\033[91mError in {result['filename']}: {result['error']}\033[0m")
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
        print(f"\033[92mTotal FASTA files checked: {len(fasta_files)}\033[0m")
        print(f"\033[92mFiles with at least one sequence match: {len(masked_files)}\033[0m")
        if masked_files:
            print("\033[92mFiles with matches:\033[0m")
            for f in sorted(masked_files):
                print(f" - {f}")
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
            print(f"\n\033[95mVerbose log saved to: {csv_path}\n")

def run_mask(args):
    os.makedirs(args.output_dir, exist_ok=True)
    sequences = load_sequences(args.sequence_file, args.output_dir)
    mask_sequences(
        args.input_dir,
        args.output_dir,
        sequences,
        verbose=args.verbose,
        dash=args.dash,
        scramble=args.scramble,
        threads=args.threads,
        dry_run=args.dry_run
    )
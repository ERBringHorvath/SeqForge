# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import os
import re
import shutil
import logging
from datetime import datetime

FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna', '.ffn', '.fas', '.faa']

def sanitize_filename(filename):
    name, ext = os.path.splitext(filename)
    name = name.replace('.', '_') #underscore
    name = re.sub(r'[-:;\s]', '_', name) #underscore
    name = re.sub(r"[()'`\"‘’“”]", '', name) #delete
    name = re.sub(r'_+', '_', name) #multi underscore to underscore
    return name + ext

def run_sanitize(args):
    logger = logging.getLogger("sanitize")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if not any(isinstance(h, logging.FileHandler) and h.baseFilename.endswith("seqforge_sanitize.log")
            for h in logger.handlers):
        file_handler = logging.FileHandler("seqforge_sanitize.log", mode='a')
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    start_time = datetime.now()
    print(f"Sanitize started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Sanitize started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    input_path = os.path.abspath(args.input)

    #arg wrangling
    if not os.path.exists(input_path):
        msg = f"Error: The specified input does not exist: {args.input}"
        print(f"\n\033[91m{msg}\033[0m")
        logger.warning(msg)
        return

    if args.in_place and args.output:
        msg = "Error: Cannot use --in-place with --sanitize-outdir. Choose one or the other"
        print(f"\n\033[91m{msg}\033[0m")
        logger.warning(msg)
        return

    if os.path.isdir(input_path) and not args.extension:
        msg = "Error: --extension (-e) is required when input is a directory."
        print(f"\n\033[91m{msg}\033[0m")
        logger.warning(msg)
        return

    if not args.in_place and not args.output:
        msg = "Error: You must specify either --in-place or --sanitize-outdir."
        print(f"\n\033[91m{msg}\033[0m")
        logger.warning(msg)
        return

    extensions = []
    if os.path.isfile(input_path):
        base_dir = os.path.dirname(input_path)
        files_to_process = [os.path.basename(input_path)]
    else:
        base_dir = input_path
        for ext in args.extension:
            if ext.lower() == 'fasta':
                extensions.extend(FASTA_EXTENSIONS)
            elif not ext.startswith('.'):
                extensions.append('.' + ext)
            else:
                extensions.append(ext)

        files_to_process = [f for f in os.listdir(base_dir) if any(f.endswith(ext) for ext in extensions)]

        if not files_to_process:
            msg = "No matching files found"
            print(f"\n\033[93m{msg}\033[0m")
            logger.warning(msg)
            return

    output_dir = os.path.abspath(args.output) if args.output else base_dir

    if not args.in_place and not args.dry_run:
        os.makedirs(output_dir, exist_ok=True)
        for f in files_to_process:
            shutil.copy2(os.path.join(base_dir, f), os.path.join(output_dir, f))

    changes_needed = []
    for f in files_to_process:
        sanitized = sanitize_filename(f)
        if sanitized != f:
            changes_needed.append((f, sanitized))

    if not changes_needed:
        msg = "No filenames require sanitization"
        print(f"\n\033[92m{msg}\033[0m")
        logger.info(msg)
        return

    print("\n\033[94mThe following filenames will be sanitized:\033[0m")
    for original, new in changes_needed:
        print(f" {original} -> {new}")
        logger.info(f"Renaming: {original} -> {new}")

    #preview of changes
    if args.dry_run:
        print("\n\033[93mDry run mode: No files were changed.\033[0m")
        logger.info("Dry run mode: No files were changed.")
        end_time = datetime.now()
        logger.info(f"Sanitization completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
        msg = f"Total runtime: {str(end_time - start_time)}"
        logger.info(msg)
        print(msg)
        return

    for original, sanitized in changes_needed:
        src_path = os.path.join(output_dir if args.output else base_dir, original)
        dst_path = os.path.join(output_dir if args.output else base_dir, sanitized)

        if args.in_place:
            os.rename(src_path, dst_path)
        else:
            os.rename(src_path, dst_path)

    msg = "Sanitization complete"
    print(f"\n\033[92m{msg}\033[0m")
    logger.info(msg)

    end_time = datetime.now()
    logger.info(f"Sanitization completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    msg = f"Total runtime: {str(end_time - start_time)}"
    logger.info(msg)
    print(msg)

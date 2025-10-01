import os
import sys
import logging
import random
import string
import hashlib
from datetime import datetime
from Bio import SeqIO

from .utils.progress import ProgressHandler

def make_suffix(seq, header, deterministic=False):
    if deterministic:
        h = hashlib.md5((seq + header).encode()).hexdigest()
        return h[:8]
    else:
        return ''.join(random.choices(string.ascii_letters + string.digits, k=5))

def run_unique_fasta_headers(args):
    logger = logging.getLogger("unique-fasta-headers")
    logger.setLevel(logging.INFO)
    if not any(isinstance(h, logging.FileHandler) and h.baseFilename.endswith("seqforge_unique-fasta-headers.log")
               for h in logger.handlers):
        fh = logging.FileHandler("seqforge_unique-fasta-headers.log")
        fh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(fh)

    start_time = datetime.now()
    print(f"Unique-Headers started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Unique-Headers started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    input_path = args.fasta_directory
    in_place = getattr(args, 'in_place', False)
    output_dir = None if in_place else getattr(args, 'output_dir', None) or getattr(args, 'output', None)

    if not in_place:
        if not output_dir:
            logger.error("No output directory provided and in-place not set")
            print("\033[91mError: Please provide --output-dir when not using --in-place")
        os.makedirs(output_dir, exist_ok=True)

    exts = {'.fa', '.fas', '.fna', '.faa', '.ffn', '.fasta'}
    fasta_files = []
    if os.path.isdir(input_path):
        for f in os.listdir(input_path):
            root, ext = os.path.splitext(f)
            if ext.lower() in exts:
                fasta_files.append(os.path.join(input_path, f))
    elif os.path.isfile(input_path):
        _, ext = os.path.splitext(input_path)
        if ext.lower() in exts:
            fasta_files = [input_path]
        else:
            logger.warning(f"Skipped {input_path}: unsupported extension.")
            print(f"\033[93mSkipped {input_path}: not a FASTA file/unsupported extension")
    else:
        logger.error(f"Input path not found: {input_path}")
        print(f"\033[91mError: Input Path not found: {input_path}")
        sys.exit(1)
    
    if not fasta_files:
        logger.warning("No FASTA files found to process")
        print("\033[91mNo FASTA files found in input path")
        return
    
    for fasta in fasta_files:
        fname = os.path.basename(fasta)
        base = os.path.splitext(fname)[0]

        try:
            records = list(SeqIO.parse(fasta, "fasta"))
        except Exception as e:
            logger.error(f"Error parsing {fasta}: {e}")
            print(f"\033[91mError parsing {fasta}: {e}")
            continue

        mode = getattr(args, 'progress', None) or 'none'
        prog = ProgressHandler(total=len(records), prefix=base, mode=mode)

        if in_place:
            dirpath = os.path.dirname(fasta) or os.getcwd()
            temp_path = os.path.join(dirpath, f"{fname}.tmp")
            target_path = temp_path
        else:
            target_path = os.path.join(output_dir, fname)
        
        try:
            with open(target_path, 'w') as out_handle:
                for rec in records:
                    # Preserve originals before mutation
                    orig_id = rec.id
                    orig_desc = rec.description  # full header (often starts with orig_id)

                    # Build a tail that excludes the leading orig_id (if present)
                    if orig_desc.startswith(orig_id):
                        tail = orig_desc[len(orig_id):].lstrip()  # drop the leading id + space(s)
                    else:
                        tail = orig_desc  # unexpected, but keep whatever was there

                    # Deterministic/random suffix using full original header
                    suffix = make_suffix(str(rec.seq), orig_desc, deterministic=args.deterministic)

                    # New stable/unique ID + preserve the original tail
                    new_id = f"{orig_id}_{base}_{suffix}"
                    rec.id = new_id
                    rec.description = tail  # keep the rest of the original header intact

                    SeqIO.write(rec, out_handle, "fasta")
                    prog.update(1, current_item=new_id)
            
            if mode in ('bar', 'verbose'):
                prog.finish()
                print()

            if in_place:
                os.replace(temp_path, fasta)
                logger.info(f"Replaced original FASTA: {fasta}")
            else:
                logger.info(f"Wrote {len(records)} records to {target_path}")

        except Exception as e:
            logger.error(f"Error writing to {target_path}: {e}")
            print(f"\033[91mError writing {target_path}: {e}\033[0m")

    end_time = datetime.now()
    print(f"Unique-Headers completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Total runtime: {str(end_time - start_time)}")
    print(f"Total runtime: {str(end_time - start_time)}")
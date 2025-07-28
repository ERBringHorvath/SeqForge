# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import os
import tempfile
import shutil
import tarfile
import zipfile

FASTA_EXTENSIONS = (
    ".fasta", ".faa", ".fna", ".ffn", ".fa", ".fas",
    ".fasta.gz", ".faa.gz", ".fna.gz", ".ffn.gz", ".fa.gz", ".fas.gz"
)

def collect_fasta_files(input_path, sanitize_flag=False):
    input_path = os.path.abspath(input_path)
    temp_dir = None

    if sanitize_flag and input_path.endswith((".zip", ".tar", ".tar.gz", ".tgz")):
        raise ValueError(
            "Sanitize (--sanitize) is incompatible with compressed archives\n"
            "Ensure all filenames inside the archive are compliant before use\n"
            "\033[92mSuggested workflow: decompress archive, run 'seqforge sanitize -i <directory> -e fasta --in-place'\033[0m\n"
            "\033[96mRe-compress with 'tar -zcvf archive.tar.gz <directory>'\033[0m"
        )
    
    def is_valid_fasta(filename):
        if filename.startswith('.') or filename.startswith('._') or filename.startswith('__MACOSX'):
            return False
        return filename.endswith(FASTA_EXTENSIONS)

    #Directory of FASTAs
    if os.path.isdir(input_path):
        fasta_files = [
            os.path.join(input_path, f)
            for f in os.listdir(input_path)
            if is_valid_fasta(f)
        ]
        if not fasta_files:
            raise ValueError(f"No FASTA files found in directory: {input_path}")
        return fasta_files, None
    
    #Single FASTA file (compressed or not)
    if os.path.isfile(input_path) and is_valid_fasta(os.path.basename(input_path)):
        return [input_path], None
    
    if os.path.isfile(input_path) and input_path.endswith((".zip", ".tar", ".tar.gz", ".tgz")):
        temp_dir = tempfile.mkdtemp(prefix="seqforge_fasta_extract_")

        try:
            if input_path.endswith(".zip"):
                with zipfile.ZipFile(input_path, 'r') as zf:
                    zf.extractall(temp_dir)
            else:
                with tarfile.open(input_path, 'r') as tf:
                    tf.extractall(temp_dir)
        except Exception as e:
            shutil.rmtree(temp_dir, ignore_errors=True)
            raise ValueError(f"Failed to extract archive {input_path}: {e}")
        
        fasta_files = []
        for root, _, files in os.walk(temp_dir):
            for file in files:
                if is_valid_fasta(file):
                    fasta_files.append(os.path.join(root, file))

        if not fasta_files:
            shutil.rmtree(temp_dir, ignore_errors=True)
            raise ValueError(f"No FASTA files found in archive: {input_path}")
        
        return fasta_files, temp_dir
    
    raise ValueError(f"No valid FASTA file(s) found at {input_path}")

def cleanup_temp_dir(temp_dir, keep=False, logger=None):
    if temp_dir and os.path.exists(temp_dir):
        if keep:
            msg = f"Temporary FASTA files retained for debugging: {temp_dir}"
            if logger:
                logger.info(msg)
            print(f"\033[93m{msg}\033[0m")
        else:
            shutil.rmtree(temp_dir, ignore_errors=True)

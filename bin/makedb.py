# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

def make_blast_db(args):

    import os, subprocess

    input_dir = args.file_directory  # Adjusted to take input directory from command line arguments
    output_dir = args.out  # Adjusted to take output directory from command line arguments

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    extensions = [".fasta", ".faa", ".fa", ".fna", ".fas"]

    for filename in os.listdir(input_dir):
        if any(filename.endswith(ext) for ext in extensions):
            input_file = os.path.join(input_dir, filename)
            base_name = os.path.splitext(filename)[0]
            output_file = os.path.join(output_dir, base_name)
            
            # Adjust dbtype based on filename extension for flexibility
            dbtype = 'nucl' if filename.endswith(('.fasta', '.fna', '.fa')) else 'prot'
            cmd = f"makeblastdb -in {input_file} -out {output_file} -dbtype {dbtype}"
            subprocess.call(cmd, shell=True)

    print("\033[92mDatabases successfully created\033[0m \n")

# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import os
import subprocess
import csv
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import islice
import re
from Bio import SeqIO
from Bio.Seq import Seq
import logging
from datetime import datetime

def execute_blast_query(data):
    (blast, db_path, query_file_path, output_file, alignment_file,
     evalue_threshold, db_name, query_file_basename, min_seq_len,
     write_alignment) = data

    # Tabular BLAST output
    cmd_tabular = (
        f"{blast} -query {query_file_path} -db {db_path} "
        f"-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen sframe' "
        f"-out {output_file} -evalue {evalue_threshold}"
    )

    if min_seq_len:
        cmd_tabular += f" -task {blast}-short -dust no"

    result = subprocess.run(cmd_tabular, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"BLAST failed: {result.stderr}")

    # Optional alignment file
    if write_alignment:
        cmd_alignment = (
            f"{blast} -query {query_file_path} -db {db_path} "
            f"-outfmt 0 -out {alignment_file} -evalue {evalue_threshold}"
        )
        subprocess.run(cmd_alignment, shell=True)

    return output_file, db_name, query_file_basename

def reverse_complement(seq):
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]

def collect_fasta_files(path, valid_exts):
    if os.path.isdir(path):
        return [os.path.join(path, f) for f in os.listdir(path) if f.endswith(tuple(valid_exts))]
    elif os.path.isfile(path) and path.endswith(tuple(valid_exts)):
        return [path]
    else:
        return []
    
def search_motif_block(rows, fasta_records, motif_regexes, logger=None):
    results = []

    for row in rows:
        gene_id = row['sseqid']
        genome = row['database']
        query = row['query_file_name']
        try:
            sstart = int(row['sstart'])
            send = int(row['send'])
            sframe = int(row['sframe'])
        except (KeyError, ValueError) as e:
            if logger:
                logger.warning(f"Skipping row with invalid numeric fields for {gene_id}: {e}")
            continue

        sequence = next(
            (seq for rec_id, seq in fasta_records.items() if rec_id.startswith(gene_id)),
            None
        )
        if not sequence:
            if logger:
                logger.warning(f"No matching sequence found for gene ID: {gene_id}")
            continue

        result_row = {
            'genome': genome,
            'query': query,
            'sseqid': gene_id,
            'sstart': sstart,
            'send': send
        }

        for idx, (motif_string, regex) in enumerate(motif_regexes, start=1):
            try:
                match = regex.search(sequence)
                result_row[f"motif_{idx}"] = match.group() if match else ''
                result_row[f"motif_{idx}_pattern"] = motif_string
            except Exception as e:
                if logger:
                    logger.error(f"Error searching motif '{motif_string}' in gene {gene_id}: {e}")
                result_row[f"motif_{idx}"] = ''
                result_row[f"motif_{idx}_pattern"] = motif_string

        # Include only if at least one match found
        if any(result_row[f"motif_{i}"] for i in range(1, len(motif_regexes)+1)):
            results.append(result_row)

    return results

def run_multiblast(args):
    logger = logging.getLogger("query")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if not any(isinstance(h, logging.FileHandler) and h.baseFilename.endswith("seqforge_query.log")
               for h in logger.handlers):
        file_handler = logging.FileHandler("seqforge_query.log", mode='a')
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    start_time = datetime.now()
    msg = f"Query stared at {start_time.strftime('%Y-%m-%d %H:%M:%S')}"
    print(msg)
    logger.info(msg)

    db_dir = args.database
    query_path = args.query_files
    fasta_path = args.fasta_directory
    threads = args.threads if args.threads else 4
    results_output_dir = args.output
    evalue_threshold = args.evalue if args.evalue is not None else 0.00001
    perc_identity_threshold = args.min_perc if args.min_perc is not None else 90
    query_coverage_threshold = args.min_cov if args.min_cov is not None else 75

    if not os.path.exists(results_output_dir):
        os.makedirs(results_output_dir)

    extensions = ['.fasta', '.fna', '.fa', '.fas', '.faa']
    db_names = set()
    tasks = []

    query_files = collect_fasta_files(query_path, extensions)
    if not query_files:
        msg = f"Error: No valid FASTA query file(s) found at {query_path}"
        print(f"\n\033[91m{msg}\033[0m")
        logger.warning(msg)

    db_basenames = set()

    if os.path.isdir(db_dir):
        for f in os.listdir(db_dir):
            if f.endswith(('.nhr', '.phr')):
                basename = os.path.splitext(f)[0]
                db_basenames.add(os.path.join(db_dir, basename))
    elif os.path.isfile(db_dir):
        if db_dir.endswith(('.nhr', '.phr')):
            db_basenames.add(os.path.splitext(db_dir)[0])
        else:
            msg = "Error: Provided file must end in .nhr or .phr"
            print(f"\033[91m{msg}\033[0m")
            logger.warning(msg)
            return
    elif os.path.exists(db_dir + ".nhr") or os.path.exists(db_dir + ".phr"):
        db_basenames.add(db_dir)
    else:
        msg = f"Error: No valid BLAST database found at {db_dir}"
        print(f"\033[91m{msg}\033[0m")
        logger.warning(msg)
        return

    db_type_map = {}
    detected_db_types = set()

    for db_path in db_basenames:
        base = os.path.basename(db_path)
        if os.path.exists(db_path + ".nhr"):
            db_type = "nucleotide"
        elif os.path.exists(db_path + ".phr"):
            db_type = "protein"
        else:
            msg = f"Warning: No valid .nhr or .phr file found for {db_path}"
            print(f"\033[93m{msg}\033[0m")
            logger.warning(msg)

        db_type_map[base] = db_type
        detected_db_types.add(db_type)

    print(f"\033[95mDetected {', '.join(detected_db_types)} database(s)\033[0m")
    logger.info(f"Detected {', '.join(detected_db_types)} database(s)")

    using_motif = False
    if args.motif:
        if "nucleotide" in detected_db_types:
            msg = "Error: --motif is only supported for protein databases (.phr)."
            print(f"\n\033[91m{msg}\033[0m")
            logger.warning(msg)
            return

        if args.nucleotide_query:
            msg = "Motif search is incompatible with --nucleotide-query. Remove one or the other"
            print(f"\n\033[91m{msg}\033[0m")
            logger.warning(msg)
            return

        if not fasta_path or not os.path.exists(fasta_path):
            msg = "Error: A valid --fasta-directory path is required with --motif"
            print(f"\n\033[91m{msg}\033[0m")
            logger.warning(msg)
            return

        valid_motifs = []
        for motif in args.motif:
            motif = motif.upper()
            if len(motif) < 4 or len(re.findall(r"[^X]", motif)) < 2:
                msg = f"Invalid motif '{motif}'. Must contain at least 4 characters and 2 defined amino acids (non-'X')"
                print(f"\n\033[91m{msg}\033[0m")
                logger.warning(msg)
                continue  # Skip this motif, but continue checking others
            valid_motifs.append(motif)

        if not valid_motifs:
            msg = "No valid motifs were provided. Exiting."
            print(f"\n\033[91m{msg}\033[0m")
            logger.warning(msg)
            return

        # Set flag and store validated motifs back into args
        using_motif = True
        args.motif = valid_motifs

    write_alignment = not args.no_alignment_files

    for query_file_path in query_files:
        query_file_basename = os.path.splitext(os.path.basename(query_file_path))[0]

        for db_file in db_basenames:
            basename = os.path.basename(db_file)
            db_type = db_type_map.get(basename)

            if not db_type:
                msg = f"Error: Could not determine type for database {basename}"
                print(f"\033[91m{msg}\033[0m")
                logger.info(msg)
                continue

            db_path = os.path.splitext(db_file)[0]
            db_name = basename

            if db_type == "nucleotide":
                blast = "blastn" if args.nucleotide_query else "tblastn"
            else:
                if args.nucleotide_query:
                    msg = f"Error: --nucleotide-query cannot be used with protein database '{db_name}'"
                    print(f"\n\033[91m{msg}\033[0m")
                    logger.warning(msg)
                    continue
                blast = "blastp"

            output_file = os.path.join(results_output_dir, f"{db_name}_{query_file_basename}_results.txt")
            alignment_file = os.path.join(results_output_dir, f"{db_name}_{query_file_basename}_alignment.txt")

            tasks.append((
                blast, db_path, query_file_path, output_file,
                alignment_file, evalue_threshold, db_name,
                query_file_basename, args.min_seq_len, write_alignment
            ))
            db_names.add((query_file_basename, db_name))

    # Run BLAST
    with ProcessPoolExecutor(max_workers=threads) as executor:
        executor.map(execute_blast_query, tasks)

    fieldnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'sframe']

    combined_results = []
    for output_file in os.listdir(results_output_dir):
        if output_file.endswith("_results.txt"):
            output_file_path = os.path.join(results_output_dir, output_file)
            with open(output_file_path, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    if len(row) == len(fieldnames):
                        task_info = next((task[6:8] for task in tasks if task[3] == output_file_path), (None, None))
                        combined_results.append(row + list(task_info))

    if not combined_results:
        print("\n \033[91mNo BLAST results found\033[0m")
        logger.info("No BLAST results found.")
        return

    df = pd.DataFrame(combined_results, columns=fieldnames + ['database', 'query_file_name'])
    df[['qstart', 'qend', 'qlen', 'pident', 'evalue', 'sframe']] = df[['qstart', 'qend', 'qlen', 'pident', 'evalue', 'sframe']].apply(pd.to_numeric, errors='coerce')
    df['query_coverage'] = ((df['qend'] - df['qstart']) / df['qlen'] * 100).round(2)

    df = df[[
        'qseqid', 'sseqid', 'database', 'query_file_name', 'pident', 'query_coverage',
        'evalue', 'bitscore', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'qlen', 'sframe']]

    df.to_csv(os.path.join(results_output_dir, "all_results.csv"), index=False)
    
    if args.motif:
        msg = f"All results saved to {results_output_dir}"
        print(f"\n\033[92m{msg}\033[0m")
        logger.info(msg)
    else:
        msg = f"All results saved to {results_output_dir}"
        print(f"\n\033[92m{msg}\033[0m\n")
        logger.info(msg)

    filtered_df = df[(df['evalue'] <= evalue_threshold) &
                     (df['pident'] >= perc_identity_threshold) &
                     (df['query_coverage'] >= query_coverage_threshold)]

    matched_pairs = set(zip(filtered_df['query_file_name'], filtered_df['database']))
    unmatched = db_names - matched_pairs
    if unmatched:
        pd.DataFrame(list(unmatched), columns=['query_file_name', 'database']).to_csv(
            os.path.join(results_output_dir, "unmatched_dbs.csv"), index=False)

    # Strongest or all filtered
    if args.report_strongest_matches:
        filtered_df.groupby(['database', 'query_file_name']).first().reset_index().to_csv(
            os.path.join(results_output_dir, "filtered_results.csv"), index=False)
    else:
        filtered_df.to_csv(os.path.join(results_output_dir, "all_filtered_results.csv"), index=False)

    if not args.keep_temp_files:
        for file in os.listdir(results_output_dir):
            if file.endswith("_results.txt"):
                try:
                    os.remove(os.path.join(results_output_dir, file))
                except Exception:
                    pass

    if using_motif:
        motif_list = args.motif
        motif_regexes = [
            (motif, re.compile(motif.replace('X', '.'), re.IGNORECASE))
            for motif in motif_list
        ]

        print(f"\n\033[95mSearching for motif pattern(s): {' '.join(motif_list)}\033[0m")
        logger.info(f"Searching for motif pattern(s): {motif_list}")
        motif_hits = []
        fasta_records = {}

        # Load FASTA records
        if os.path.isdir(fasta_path):
            for file in os.listdir(fasta_path):
                if any(file.endswith(ext) for ext in extensions):
                    for record in SeqIO.parse(os.path.join(fasta_path, file), "fasta"):
                        fasta_records[record.id] = str(record.seq).upper()
        elif os.path.isfile(fasta_path):
            for record in SeqIO.parse(fasta_path, "fasta"):
                fasta_records[record.id] = str(record.seq).upper()

        # Chunk generator
        from itertools import islice
        def chunkify(iterable, n):
            it = iter(iterable)
            return iter(lambda: list(islice(it, n)), [])

        rows = [row for _, row in df.iterrows()]
        row_chunks = chunkify(rows, max(1, len(rows) // threads))

        # Run parallel motif search
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            futures = [
                executor.submit(search_motif_block, chunk, fasta_records, motif_regexes)
                for chunk in row_chunks
            ]
            for future in futures:
                result = future.result()
                if result:
                    motif_hits.extend(result)

        if motif_hits:
            motif_df = pd.DataFrame(motif_hits)
            base_cols = ['genome', 'query', 'sseqid', 'sstart', 'send']
            motif_cols = []
            for idx in range(1, len(motif_regexes) + 1):
                motif_cols.extend([f"motif_{idx}", f"motif_{idx}_pattern"])

            ordered_cols = base_cols + [col for col in motif_cols if col in motif_df.columns]
            motif_df = motif_df[ordered_cols]

            for idx in range(1, len(motif_regexes) + 1):
                count = motif_df[f"motif_{idx}"].astype(bool).sum()
                msg = f"Motif {idx} ({motif_regexes[idx-1][0]}): {count} hits"
                print(f"\033[92m{msg}\033[0m")
                logger.info(msg)
            
            motif_df.to_csv(os.path.join(results_output_dir, "motif_matches.csv"), index=False)
            print("\n\033[92mMotif matches saved\033[0m\n")
            logger.info("Motif matches saved")
        else:
             print("\n\033[91mNo motif matches found\033[0m")
             logger.warning("No motif matches found")

    if args.visualize:
        from visualize import run_visualization
        run_visualization(filtered_df, args)

    if args.visualize and args.motif:
        if motif_hits:
            from visualize import run_sequence_logo
            motif_df = pd.DataFrame(motif_hits)
            run_sequence_logo(motif_df, args)

    end_time = datetime.now()
    print(f"\nQuery completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Query completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    msg = f"Total runtime: {str(end_time - start_time)}"
    logger.info(msg)
    print(msg)            

    return filtered_df

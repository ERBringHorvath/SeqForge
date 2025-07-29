# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import os
import subprocess
import csv
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from itertools import islice
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from datetime import datetime
import logging

from utils.file_handler import collect_fasta_files, cleanup_temp_dir

FASTA_EXTS = ('.fa', '.fna', '.fas', '.ffn', '.fasta', '.faa')

def chunkify(iterable, n):
    #Split iterable into chunks of size n for motif mining
    it = iter(iterable)
    return iter(lambda: list(islice(it, n)), [])

def export_motif_fasta(motif_df, output_dir, fasta_path=None, motif_only=False):
    if motif_df.empty:
        print("\033[93mNo motif matches to write FASTA files for.\033[0m")
        return

    fasta_records = {}
    if not motif_only and fasta_path:
        try:
            fasta_files, temp_dir = collect_fasta_files(fasta_path)
        except ValueError:
            fasta_files, temp_dir = [], None

        for file in fasta_files:
            for record in SeqIO.parse(file, "fasta"):
                fasta_records[record.id] = str(record.seq).upper()

        cleanup_temp_dir(temp_dir, keep=False)

    motif_cols = [col for col in motif_df.columns if col.startswith("motif_") and not col.endswith("_pattern")]

    for idx, motif_col in enumerate(motif_cols, start=1):
        pattern_col = f"motif_{idx}_pattern"
        motif_rows = motif_df[motif_col].dropna()
        motif_rows = motif_rows[motif_rows.str.strip() != ""]

        if motif_rows.empty:
            continue

        motif_pattern = motif_df[pattern_col].dropna().iloc[0]
        fasta_out = os.path.join(output_dir, f"{motif_pattern}_matches.faa")

        with open(fasta_out, "w") as f:
            for i, row in enumerate(motif_df.loc[motif_rows.index].itertuples(), start=1):
                header = f">{row.query}_{row.sstart}_{row.send}_hit{i}"

                if motif_only:
                    sequence = getattr(row, motif_col)
                else:
                    #Extract full BLAST-aligned region from the FASTA record
                    seq = fasta_records.get(row.sseqid)
                    if not seq:
                        continue
                    low, high = sorted([row.sstart, row.send])
                    sequence = seq[low - 1:high]  #1-based to 0-based adjustment

                f.write(f"{header}\n{sequence}\n")

        print(f"\033[92mWrote {len(motif_rows)} matches to {fasta_out}\033[0m")

def execute_blast_query(data):
    (blast, db_path, query_file_path, output_file, alignment_file,
     evalue_threshold, db_name, query_file_basename, min_seq_len,
     write_alignment) = data

    cmd_tabular = (
        f"{blast} -query {query_file_path} -db {db_path} "
        f"-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
        f"evalue bitscore qlen sframe' "
        f"-out {output_file} -evalue {evalue_threshold}"
    )
    if min_seq_len:
        cmd_tabular += f" -task {blast}-short -dust no"

    result = subprocess.run(cmd_tabular, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"\033[91mBLAST failed:\033[0m {result.stderr}")

    if write_alignment:
        cmd_alignment = (
            f"{blast} -query {query_file_path} -db {db_path} "
            f"-outfmt 0 -out {alignment_file} -evalue {evalue_threshold}"
        )
        subprocess.run(cmd_alignment, shell=True)

    return output_file, db_name, query_file_basename

def search_motif_block(rows, fasta_records, motif_regexes):
    import pandas as pd
    import os

    results, warnings = [], []
    matched_restrictions = {idx: False for idx, (_, _, target) in enumerate(motif_regexes, start=1) if target}

    for row in rows:
        gene_id = row['sseqid']
        genome = row['database']
        query_file = row['query_file_name']
        query_basename = os.path.splitext(query_file)[0]

        try:
            sstart = int(row['sstart'])
            send = int(row['send'])
        except (KeyError, ValueError) as e:
            warnings.append(f"Skipping row with invalid numeric fields for {gene_id}: {e}")
            continue

        sequence = next((seq for rec_id, seq in fasta_records.items() if rec_id.startswith(gene_id)), None)
        if not sequence:
            warnings.append(f"No matching sequence found for gene ID: {gene_id}")
            continue

        result_row = {
            'genome': genome,
            'query': query_file,
            'sseqid': gene_id,
            'sstart': sstart,
            'send': send
        }

        #Search motifs
        for idx, (motif_string, regex, target_basename) in enumerate(motif_regexes, start=1):
            if target_basename and query_basename != target_basename:
                result_row[f"motif_{idx}"] = ''
                result_row[f"motif_{idx}_pattern"] = motif_string
                continue

            try:
                matches = [m.group() for m in regex.finditer(sequence)]
                result_row[f"motif_{idx}"] = ",".join(matches) if matches else ''
                result_row[f"motif_{idx}_pattern"] = motif_string
                if matches and target_basename:
                    matched_restrictions[idx] = True
            except Exception as e:
                warnings.append(f"Error searching motif '{motif_string}' in gene {gene_id}: {e}")
                result_row[f"motif_{idx}"] = ''
                result_row[f"motif_{idx}_pattern"] = motif_string

        #Only keep rows with at least one motif hit
        if any(result_row[f"motif_{i}"] for i in range(1, len(motif_regexes) + 1)):
            results.append(result_row)

    if results:
        df = pd.DataFrame(results)
        final_groups = []

        for gene_id, group in df.groupby('sseqid'):
            group = group.copy()
            motif_cols = [c for c in group.columns if c.startswith('motif_') and not c.endswith('_pattern')]

            for col in motif_cols:
                #Skip columns with no commas entirely
                if not any(',' in str(x) for x in group[col]):
                    continue

                split_motifs = [x.split(',') if x else [] for x in group[col]]
                cleaned = []
                for row_idx, motifs in enumerate(split_motifs):
                    #Only adjust if this row had a comma-separated list
                    if len(motifs) > 1:
                        cleaned.append(motifs[row_idx] if row_idx < len(motifs) else '')
                    else:
                        cleaned.append(motifs[0] if motifs else '')
                group[col] = cleaned

            final_groups.append(group)

        df = pd.concat(final_groups, ignore_index=True)
        results = df.to_dict(orient='records')

    return results, warnings

def run_motif_search(df, fasta_path, motif_regexes, results_output_dir, threads):
    fasta_records = {}
    try:
        fasta_files, temp_dir = collect_fasta_files(fasta_path)
    except ValueError as e:
        print(f"\n\033[91mError: {e}\033[0m")
        return pd.DataFrame()
    
    for file in fasta_files:
        for record in SeqIO.parse(file, "fasta"):
            fasta_records[record.id] = str(record.seq).upper()
    
    cleanup_temp_dir(temp_dir, keep=False)

    #Parallel search
    rows = [row for _, row in df.iterrows()]
    row_chunks = chunkify(rows, max(1, len(rows) // max(1, threads)))

    results, warnings = [], []
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(search_motif_block, chunk, fasta_records, motif_regexes)
                   for chunk in row_chunks]
        for future in futures:
            res, warns = future.result()
            results.extend(res)
            warnings.extend(warns)

    #Log warnings once
    for w in warnings:
        print(f"\033[93m{w}\033[0m")

    if results:
        motif_df = pd.DataFrame(results)
        base_cols = ['genome', 'query', 'sseqid', 'sstart', 'send']

        #Ensure motif columns always exist
        for idx in range(1, len(motif_regexes) + 1):
            for col in (f"motif_{idx}", f"motif_{idx}_pattern"):
                if col not in motif_df.columns:
                    motif_df[col] = ''

        motif_cols = [c for idx in range(1, len(motif_regexes)+1)
                      for c in (f"motif_{idx}", f"motif_{idx}_pattern")]
        motif_df = motif_df[base_cols + motif_cols]

        motif_df.to_csv(os.path.join(results_output_dir, "motif_matches.csv"), index=False)
        print(f"\n\033[92mMotif matches saved to motif_matches.csv\033[0m")
    else:
        print("\n\033[91mNo motif matches found\033[0m")
        motif_df = pd.DataFrame()
    
    return motif_df

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
    print(f"Query started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Query started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    db_dir = args.database
    query_path = args.query_files
    fasta_path = args.fasta_directory
    results_output_dir = args.output
    threads = args.threads if args.threads else 4

    evalue_threshold = args.evalue if args.evalue is not None else 1e-5
    perc_identity_threshold = args.min_perc if args.min_perc is not None else 90
    query_coverage_threshold = args.min_cov if args.min_cov is not None else 75

    if not os.path.exists(results_output_dir):
        os.makedirs(results_output_dir)

    try:
        query_files, temp_dir = collect_fasta_files(query_path)
    except ValueError as e:
        print(f"\n\033[91mError: {e}\033[0m")
        logger.error(str(e))
        return
    
    if not query_files:
        print(f"\033[91mError: No valid FASTA file(s) found at {query_path}\033[0m")
        logger.warning("No valid FASTA files found")
        return

    db_basenames = set()
    if os.path.isdir(db_dir):
        for f in os.listdir(db_dir):
            if f.endswith(('.nhr', '.phr')):
                db_basenames.add(os.path.splitext(os.path.join(db_dir, f))[0])
    elif os.path.isfile(db_dir):
        if db_dir.endswith(('.nhr', '.phr')):
            db_basenames.add(os.path.splitext(db_dir)[0])
        else:
            print(f"\033[91mError: Provided database file must end in .nhr or .phr\033[0m")
            return
    elif os.path.exists(db_dir + ".nhr") or os.path.exists(db_dir + ".phr"):
        db_basenames.add(db_dir)
    else:
        print(f"\033[91mError: No valid BLAST database found at {db_dir}\033[0m")
        logger.warning(f"No valid BLAST database found at {db_dir}")
        return

    #Determine DB types
    db_type_map, detected_db_types = {}, set()
    for db_path in db_basenames:
        base = os.path.basename(db_path)
        if os.path.exists(db_path + ".nhr"):
            db_type = "nucleotide"
        elif os.path.exists(db_path + ".phr"):
            db_type = "protein"
        else:
            print(f"\033[93mWarning: No valid .nhr or .phr for {db_path}\033[0m")
            logger.warning(f"No valid DB index for {db_path}")
            continue
        db_type_map[base] = db_type
        detected_db_types.add(db_type)

    print(f"\033[95mDetected {', '.join(detected_db_types)} database(s)\033[0m")
    logger.info(f"Detected {', '.join(detected_db_types)} database(s)")

    #Validate and preprocess motifs
    motif_regexes, using_motif = [], False
    if args.motif:
        if "nucleotide" in detected_db_types:
            print("\033[91mError: --motif only supported for protein databases (.phr)\033[0m")
            return
        if args.nucleotide_query:
            print("\033[91mError: --motif incompatible with --nucleotide-query\033[0m")
            return
        if not fasta_path or not os.path.exists(fasta_path):
            print("\033[91mError: A valid --fasta-directory path is required with --motif\033[0m")
            return

        #Prepare motif list, allow {basename} filtering
        query_basenames = [os.path.splitext(os.path.basename(f))[0] for f in query_files]
        valid_motifs = []

        for raw_motif in args.motif:
            target_basename = None
            if (('{' in raw_motif and raw_motif.endswith(('}', ']'))) or 
                ('[' in raw_motif and raw_motif.endswith((']', '}')))):
                try:
                    if '{' in raw_motif:
                        motif_str, target_basename = raw_motif.split('{', 1)
                        target_basename = target_basename.rstrip('}]')
                    else:
                        motif_str, target_basename = raw_motif.split('[', 1)
                        target_basename = target_basename.rstrip('}]')

                    motif = motif_str.upper()

                    if not target_basename or target_basename not in query_basenames:
                        print(f"\033[93mWarning: Motif '{motif}' restricted to '{target_basename}', "
                              f"but no query file matched. Skipping\033[0m")
                        continue

                except ValueError:
                    print(f"\033[91mInvalid motif format '{raw_motif}' "
                          f"use MOTIF{{basename}} or MOTIF[basename]\033[0m")
                    continue

            else:
                motif = raw_motif.upper()
                target_basename = None

            if len(motif) < 4 or len(re.findall(r"[^X]", motif)) < 2:
                print(f"\033[91mInvalid motif '{motif}': must have ≥ 4 characters and ≥ 2 non-'X'\033[0m")
                continue

            valid_motifs.append((motif, re.compile(motif.replace("X", "."), re.IGNORECASE), target_basename))

        if not valid_motifs:
            print("\033[91mNo valid motifs provided. Exiting.\033[0m")
            return

        motif_regexes = valid_motifs
        using_motif = True

    #Build BLAST tasks
    db_names, tasks = set(), []
    write_alignment = not args.no_alignment_files
    for query_file_path in query_files:
        query_file_basename = os.path.splitext(os.path.basename(query_file_path))[0]
        for db_file in db_basenames:
            basename = os.path.basename(db_file)
            db_type = db_type_map.get(basename)
            if not db_type:
                print(f"\033[91mError: Could not determine DB type for {basename}\033[0m")
                continue

            blast = "blastn" if db_type == "nucleotide" and args.nucleotide_query else (
                    "tblastn" if db_type == "nucleotide" else "blastp")
            if args.nucleotide_query and db_type == "protein":
                print(f"\033[91mError: --nucleotide-query not valid for protein DB '{db_file}'\033[0m")
                continue

            output_file = os.path.join(results_output_dir, f"{basename}_{query_file_basename}_results.txt")
            alignment_file = os.path.join(results_output_dir, f"{basename}_{query_file_basename}_alignment.txt")

            tasks.append((blast, db_file, query_file_path, output_file, alignment_file,
                          evalue_threshold, basename, query_file_basename, args.min_seq_len,
                          write_alignment))
            db_names.add((query_file_basename, basename))

    with ProcessPoolExecutor(max_workers=threads) as executor:
        executor.map(execute_blast_query, tasks)

    fieldnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'sframe']
    combined_results = []
    for output_file in os.listdir(results_output_dir):
        if output_file.endswith("_results.txt"):
            path = os.path.join(results_output_dir, output_file)
            with open(path, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    if len(row) == len(fieldnames):
                        task_info = next((task[6:8] for task in tasks if task[3] == path), (None, None))
                        combined_results.append(row + list(task_info))

    if not combined_results:
        print("\033[91mNo BLAST results found.\033[0m")
        return

    df = pd.DataFrame(combined_results, columns=fieldnames + ['database', 'query_file_name'])
    df[['qstart', 'qend', 'qlen', 'pident', 'evalue', 'sframe']] = (
        df[['qstart', 'qend', 'qlen', 'pident', 'evalue', 'sframe']].apply(pd.to_numeric, errors='coerce')
    )
    df['query_coverage'] = ((df['qend'] - df['qstart']) / df['qlen'] * 100).round(2)

    df = df[['qseqid', 'sseqid', 'database', 'query_file_name', 'pident', 'query_coverage',
             'evalue', 'bitscore', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
             'sstart', 'send', 'qlen', 'sframe']]
    df.to_csv(os.path.join(results_output_dir, "all_results.csv"), index=False)
    print(f"\033[92mSaved all_results.csv\033[0m")

    filtered_df = df[(df['evalue'] <= evalue_threshold) &
                     (df['pident'] >= perc_identity_threshold) &
                     (df['query_coverage'] >= query_coverage_threshold)]

    if args.report_strongest_matches:
        filtered_df.groupby(['database', 'query_file_name']).first().reset_index().to_csv(
            os.path.join(results_output_dir, "filtered_results.csv"), index=False)
    else:
        filtered_df.to_csv(os.path.join(results_output_dir, "all_filtered_results.csv"), index=False)

    if not args.keep_temp_files:
        for file in os.listdir(results_output_dir):
            if file.endswith("_results.txt"):
                os.remove(os.path.join(results_output_dir, file))

    motif_df = pd.DataFrame()
    if using_motif:
        print(f"\033[95mSearching for motifs: {' '.join(m[0] for m in motif_regexes)}\033[0m")
        motif_df = run_motif_search(df, fasta_path, motif_regexes, results_output_dir, threads)

        #Handle restricted-target warnings only if hits exist
        if not motif_df.empty:
            restricted_targets = {target for _, _, target in motif_regexes if target}
            seen_queries = set(motif_df['query'].unique())
            for target in restricted_targets:
                if target not in seen_queries:
                    msg = f"Warning: Motif restricted to '{target}' had no hits in any query results."
                    print(f"\033[93m{msg}\033[0m")
                    logger.warning(msg)

    if using_motif and not motif_df.empty and args.motif_fasta_out:
        export_motif_fasta(motif_df, results_output_dir, 
                           fasta_path=fasta_path,
                           motif_only=args.motif_only)

    if args.visualize:
        from visualize import run_visualization
        run_visualization(filtered_df, args)
        if using_motif and not motif_df.empty:
            from visualize import run_sequence_logo
            run_sequence_logo(motif_df, args)

    cleanup_temp_dir(temp_dir, keep=getattr(args, 'keep_temp_files', False), logger=logger)

    end_time = datetime.now()
    print(f"Query completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Total runtime: {str(end_time - start_time)}")
    print(f"Total runtime: {str(end_time - start_time)}")
    return filtered_df

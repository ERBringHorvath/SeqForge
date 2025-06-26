import os
import subprocess
import csv
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

def execute_blast_query(data):
    (blast, db_path, query_file_path, output_file, alignment_file,
     evalue_threshold, db_name, query_file_base_name, min_seq_len, write_alignment) = data

    # Tabular output
    cmd_tabular = (
        f"{blast} -query {query_file_path} -db {db_path} "
        f"-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' "
        f"-out {output_file} -evalue {evalue_threshold}"
    )

    if min_seq_len:
        cmd_tabular += f" -task {blast}-short -dust no"

    subprocess.run(cmd_tabular, shell=True)

    # Optional alignment output
    if write_alignment:
        cmd_alignment = (
            f"{blast} -query {query_file_path} -db {db_path} "
            f"-outfmt 0 -out {alignment_file} -evalue {evalue_threshold}"
        )
        subprocess.run(cmd_alignment, shell=True)

    return output_file, db_name, query_file_base_name

def run_multiblast(args):
    blast = 'blastn' if args.nucleotide_query else 'tblastn'
    db_dir = args.database
    query_path = args.query_files
    threads = args.threads if args.threads else 1
    results_output_dir = args.output
    evalue_threshold = args.evalue if args.evalue is not None else 0.00001
    perc_identity_threshold = args.perc if args.perc is not None else 90
    query_coverage_threshold = args.cov if args.cov is not None else 75

    # Handle directories
    if not os.path.exists(db_dir):
        print(f"\n \033[91mNo BLAST databases found at {db_dir}\033[0m")
        return

    if not os.path.exists(query_path):
        print(f"No query files found at {query_path}")
        return

    if not os.path.exists(results_output_dir):
        os.makedirs(results_output_dir)

    extensions = ['.fasta', '.fna', '.fa', '.fas', '.faa']
    fieldnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen']

    write_alignment = not args.no_alignment_files

    tasks = []
    db_names = set()
    for query_file in os.listdir(query_path):
        if any(query_file.endswith(ext) for ext in extensions):
            query_file_path = os.path.join(query_path, query_file)
            query_file_base_name = os.path.splitext(query_file)[0]

            for file_name in os.listdir(db_dir):
                if file_name.endswith(".nhr"):
                    db_name = os.path.splitext(file_name)[0]
                    db_path = os.path.join(db_dir, db_name)
                    output_file_path = os.path.join(results_output_dir, f"{db_name}_{query_file_base_name}_results.txt")
                    alignment_file_path = os.path.join(results_output_dir, f"{db_name}_{query_file_base_name}_alignment.txt")
                    tasks.append((
                        blast, db_path, query_file_path, output_file_path,
                        alignment_file_path, evalue_threshold, db_name,
                        query_file_base_name, args.min_seq_len, write_alignment
                    ))
                    db_names.add((query_file_base_name, db_name))

    # Run BLAST queries in parallel
    with ProcessPoolExecutor(max_workers=threads) as executor:
        executor.map(execute_blast_query, tasks)

    # Combine results
    combined_results = []
    for output_file in os.listdir(results_output_dir):
        if output_file.endswith("_results.txt"):
            output_file_path = os.path.join(results_output_dir, output_file)
            with open(output_file_path, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    if len(row) == len(fieldnames):
                        task_info = next((task[5:7] for task in tasks if task[3] == output_file_path), (None, None))
                        combined_results.append(row + list(task_info))

    if not combined_results:
        print("\n \033[91mNo BLAST results found. Please check input files and parameters\033[0m")
        return

    df = pd.DataFrame(combined_results, columns=fieldnames + ['database', 'query_file_name'])
    df[['qstart', 'qend', 'qlen', 'pident', 'evalue']] = df[['qstart', 'qend', 'qlen', 'pident', 'evalue']].apply(pd.to_numeric, errors='coerce')
    df['query_coverage'] = ((df['qend'] - df['qstart']) / df['qlen'] * 100).round(2)

    column_order = [
        'qseqid', 'sseqid', 'database', 'query_file_name', 'pident', 'query_coverage', 'evalue', 'bitscore',
        'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'qlen'
    ]

    df = df[column_order]

    output_csv_path = os.path.join(results_output_dir, "all_results.csv")
    df.to_csv(output_csv_path, index=False)
    print(f"\n \033[92mAll results saved in {output_csv_path}\033[0m")

    filtered_df = df[(df['evalue'] <= evalue_threshold) & (df['pident'] >= perc_identity_threshold) & (df['query_coverage'] >= query_coverage_threshold)]

   # Identify databases with no hits
    matched_db_pairs = set(zip(filtered_df['query_file_name'], filtered_df['database']))
    unmatched_db_pairs = db_names - matched_db_pairs

    if unmatched_db_pairs:
        unmatched_df = pd.DataFrame(list(unmatched_db_pairs), columns=['query_file_name', 'database'])
        unmatched_output_path = os.path.join(results_output_dir, "unmatched_dbs.csv")
        unmatched_df.to_csv(unmatched_output_path, index=False)
        print(f"\n \033[93mUnmatched databases saved in {unmatched_output_path}\033[0m")

    if args.report_strongest_matches:
        strongest_hits = filtered_df.groupby(['database', 'query_file_name']).first().reset_index()
        output_csv_path = os.path.join(results_output_dir, "filtered_results.csv")
        strongest_hits.to_csv(output_csv_path, index=False)
        print(f"\n \033[92mStrongest matches stored in {output_csv_path}\033[0m")
    else:
        output_csv_path = os.path.join(results_output_dir, "all_filtered_results.csv")
        filtered_df.to_csv(output_csv_path, index=False)
        print(f"\n \033[92mFiltered results stored in {output_csv_path}\033[0m")

    if not args.keep_temp_files:
        deleted_count = 0
        print("\n \033[94mDeleting temporary files....\033[0m")
        for file in os.listdir(results_output_dir):
            if file.endswith("_results.txt"):
                try:
                    os.remove(os.path.join(results_output_dir, file))
                    deleted_count += 1
                except Exception as e:
                    print(f"\n \033[91mWarning: Failed to delete {file} - {e}\033[0m")
        print(f" \033[92mDeleted {deleted_count} temporary file(s)\033[0m")

    return filtered_df

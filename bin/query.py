import os
import subprocess
import csv
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import re
from Bio import SeqIO
from Bio.Seq import Seq

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

def run_multiblast(args):
    db_dir = args.database
    query_path = args.query_files
    fasta_path = args.fasta_directory
    threads = args.threads if args.threads else 4
    results_output_dir = args.output
    evalue_threshold = args.min_evalue if args.min_evalue is not None else 0.00001
    perc_identity_threshold = args.min_perc if args.min_perc is not None else 90
    query_coverage_threshold = args.min_cov if args.min_cov is not None else 75

    if not os.path.exists(results_output_dir):
        os.makedirs(results_output_dir)

    extensions = ['.fasta', '.fna', '.fa', '.fas', '.faa']
    db_names = set()
    tasks = []

    # Parse queries
    query_files = []
    if os.path.isdir(query_path):
        for file in os.listdir(query_path):
            if any(file.endswith(ext) for ext in extensions):
                query_files.append(os.path.join(query_path, file))
    elif os.path.isfile(query_path) and any(query_path.endswith(ext) for ext in extensions):
        query_files.append(query_path)
    else:
        print(f"\n\033[91mError: No valid FASTA query file(s) found at {query_path}\033[0m")
        return

    db_type_map = {}
    detected_db_types = set()

    for file_name in os.listdir(db_dir):
        basename, ext = os.path.splitext(file_name)
        if ext in [".nhr", ".phr"]:
            db_type = "nucleotide" if ext == ".nhr" else "protein"
            detected_db_types.add(db_type)
    
    print(f"\033[95mDetected {db_type} database(s)\033[0m")

    # Validate motif input
    using_motif = False
    if args.motif:
        if "nucleotide" in detected_db_types:
            print("\n\033[91mError: --motif is only supported for protein databases (.phr).\033[0m")
            return
        motif = args.motif.upper()
        if len(motif) < 4 or len(re.findall(r"[^X]", motif)) < 2:
            print("\n\033[91mInvalid motif. Must contain at least 4 characters and 2 defined amino acids (non-'X')\033[0m")
            return
        if args.nucleotide_query:
            print("\n\033[91mMotif search is incompatible with --nucleotide-query. Remove one or the other\033[0m")
            return
        if not fasta_path or not os.path.exists(fasta_path):
            print("\n\033[91mError: A valid --fasta-directory path is required with --motif\033[0m")
            return
        using_motif = True
        motif_regex = motif.replace("X", ".")

    write_alignment = not args.no_alignment_files

    # Prepare BLAST tasks
    for query_file_path in query_files:
        query_file_basename = os.path.splitext(os.path.basename(query_file_path))[0]

        for file_name in os.listdir(db_dir):
            basename, ext = os.path.splitext(file_name)
            if ext not in [".nhr", ".phr"]:
                continue

            db_type = "nucleotide" if ext == ".nhr" else "protein"
            db_type_map[basename] = db_type

            db_path = os.path.join(db_dir, basename)
            db_name = basename

            if db_type == "nucleotide":
                blast = "blastn" if args.nucleotide_query else "tblastn"
            else:
                if args.nucleotide_query:
                    print(f"\n\033[91mError: --nucleotide-query cannot be used with protein database '{db_name}'\033[0m")
                    continue
                blast = 'blastp'

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

    # Parse results
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
        print("\n \033[91mNo BLAST results found.\033[0m")
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
        print(f"\n\033[92mAll results saved to {results_output_dir}.\033[0m")
    else:
        print(f"\n\033[92mAll results saved to {results_output_dir}.\033[0m\n")

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

    # Motif search logic
    if using_motif:
        print(f"\n\033[95mSearching for motif pattern: {motif}\033[0m")
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

        for _, row in filtered_df.iterrows():
            gene_id = row['sseqid']
            genome = row['database']
            query = row['query_file_name']
            sstart, send = int(row['sstart']), int(row['send'])
            sframe = int(row['sframe'])

            db_type = db_type_map.get(genome)
            sequence = next((seq for rec_id, seq in fasta_records.items() if rec_id.startswith(gene_id)), None)
            if not sequence:
                continue

            try:
                match = re.search(motif_regex, sequence)
                if match:
                    motif_hits.append({
                        'genome': genome,
                        'query': query,
                        'motif_match': match.group()
                    })
                else:
                    print("\n\033[91mError parsing amino acid motif\033[0m")

            # try:
            #     if db_type == 'nucleotide':
            #         start, end = min(sstart, send) - 1, max(sstart, send)
            #         if start < 0 or end > len(sequence):
            #             print(f"Skipping out-of-bounds: {start}-{end} on seq len {len(sequence)} ({gene_id})")
            #             continue

            #         nt_seq = sequence[start:end]
            #         if sframe < 0:
            #             nt_seq = str(Seq(nt_seq).reverse_complement())

            #         offset = abs(sframe) - 1
            #         nt_seq = nt_seq[offset:]

            #         pad_len = len(nt_seq) % 3
            #         if pad_len != 0:
            #             nt_seq += "N" * (3 - pad_len)

            #         aa_seq = str(Seq(nt_seq).translate(to_stop=False))
            #         match = re.search(motif_regex, aa_seq)

            #         if match:
            #             motif_hits.append({
            #                 'genome': genome,
            #                 'query': query,
            #                 'motif_match': match.group()
            #             })

            #     else:  # blastp
            #         match = re.search(motif_regex, sequence)
            #         if match:
            #             motif_hits.append({
            #                 'genome': genome,
            #                 'query': query,
            #                 'motif_match': match.group()
            #             })

            except Exception as e:
                print(f"Error matching {gene_id}: {e}")

        if motif_hits:
            pd.DataFrame(motif_hits).to_csv(
                os.path.join(results_output_dir, "motif_matches.csv"), index=False)
            print(f"\n\033[92mMotif matches saved.\033[0m")
        else:
            print("\n\033[91mNo motif matches found.\033[0m")

    if not args.keep_temp_files:
        for file in os.listdir(results_output_dir):
            if file.endswith("_results.txt"):
                try:
                    os.remove(os.path.join(results_output_dir, file))
                except Exception:
                    pass

    return filtered_df
import os
import subprocess
import csv
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

##Define blast command chain
def execute_blast_query(data):
        blast, db_path, query_file_path, output_file, e_value_threshold = data  # Adjusted parameters
        cmd = f"{blast} -query {query_file_path} -db {db_path} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -out {output_file} -evalue {e_value_threshold}"
        subprocess.run(cmd, shell=True)
        print(f"Completed BLAST for {query_file_path}")
        return output_file  # Return individual output file path for processing

##Blast and process results
def run_multiblastp(args):
    """
    Runs BLAST queries in parallel
    """

    ##Define arguments to be passed from user
    blast = args.method
    db_dir = args.database
    query_path = args.query_files
    threads = args.threads
    results_output_dir = args.output
    e_value_threshold = args.evalue

    ##Handle directories
    if not os.path.exists(db_dir):
        print(f"No BLAST databases found at {db_dir}")
        return
    
    if not os.path.exists(query_path):
        print(f"No query files found at {query_path}")
        return
    
    if not os.path.exists(results_output_dir):
        os.makedirs(results_output_dir)

    ##Create function to process output files
    def process_output_file(output_file, db_name, query_file_base_name, e_value_threshold, fieldnames):
        results = []
        with open(output_file, 'r') as f:
            reader = csv.DictReader(f, fieldnames=fieldnames[:-3], delimiter='\t')
            for row in reader:
                if float(row['evalue']) <= e_value_threshold:
                    row['database'] = db_name
                    row['query_file_name'] = query_file_base_name
                    row['file_path'] = output_file
                    results.append(row)  # Append row to results list
        return results

    extensions = ['.fasta', '.fna', '.fa', '.fas', '.faa']
    fieldnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'file_path', 'database', 'query_file_name']

    ##Define tasks
    tasks = []
    for query_file in os.listdir(query_path):
        if any(query_file.endswith(ext) for ext in extensions):
            query_file_path = os.path.join(query_path, query_file)
            query_file_base_name = os.path.splitext(query_file)[0]

            for file_name in os.listdir(db_dir):
                if file_name.endswith(".nhr"):
                    db_name = os.path.splitext(file_name)[0]
                    db_path = os.path.join(db_dir, db_name)
                    output_file_path = os.path.join(results_output_dir, f"{db_name}_{query_file_base_name}_results.txt")  # Corrected variable name
                    tasks.append((blast, db_path, query_file_path, output_file_path, e_value_threshold))

    #Parallel execution call
    with ProcessPoolExecutor(max_workers=threads) as executor:
        output_files = list(executor.map(execute_blast_query, tasks))

    ##Combine results
    results = []
    for output_file in output_files:
        processed_results = process_output_file(output_file, db_name, query_file_base_name, e_value_threshold, fieldnames)
        results.extend(processed_results)

    df = pd.DataFrame(results)

    ##Rename headers
    column_names_mapping = {
    'database': 'database',
    'query_file_name': 'query_file_name',
    'evalue': 'evalue',
    'bitscore': 'bitscore',
    'pident': 'percent_identity',
    'sseqid': 'location',
    'sstart': 'subject_start',
    'send': 'subject_end',
    'qseqid': 'query_ID',
    'length': 'alignment_length',
    'mismatch': 'mismatches',
    'gapopen': 'gap_opens',
    'qstart': 'query_start',
    'qend': 'query_end',
    'file_path': 'result_file_path'
}

    #Define column order
    ordered_fieldnames = [
        'database',        
        'query_file_name',      
        'evalue',              
        'bitscore',
        'percent_identity',            
        'location',          
        'subject_start',        
        'subject_end',         
        'query_ID',
        'alignment_length',    
        'mismatches',          
        'gap_opens',
        'query_start',
        'query_end',
        'result_file_path'
    ]
    
    df = df.rename(columns=column_names_mapping)
    df = df[ordered_fieldnames]

    # Save DataFrame to CSV
    df.to_csv(os.path.join(results_output_dir, "results.csv"), index=False)

    print(f"\n \033[92mQueries Complete and stored in {results_output_dir}\033[0m \n")

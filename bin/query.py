def run_blast_query(args):
    """
    Runs BLAST queries against specified databases and compiles the results.
    """

    import os, subprocess, csv
    import pandas as pd

    blast = args.method
    db_dir = args.database
    query_path = args.query_files
    results_output = args.output
    e_value_threshold = args.evalue if args.evalue else 0.001
    report_only_lowest_evalue = args.report_only_lowest_evalue

    if not os.path.exists(db_dir):
        print(f"No BLAST databases found at {db_dir}")
        return

    if not os.path.exists(query_path):
        print(f"No query files found at {query_path}")
        return

    if not os.path.exists(results_output):
        os.makedirs(results_output)

    extensions = ['.fasta', '.fna', '.fa', '.fas', '.faa']
    fieldnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'file_path', 'database', 'query_file_name']
    results = []

    print("\n \033[92mProcessing queries...\033[0m")

    print(f"\n \033[93mReport only top hit: {report_only_lowest_evalue}\033[0m")

    for query_file in os.listdir(query_path):
        if any(query_file.endswith(ext) for ext in extensions):
            query_file_path = os.path.join(query_path, query_file)
            query_file_base_name = os.path.splitext(query_file)[0]

            for file_name in os.listdir(db_dir):
                if file_name.endswith(".nhr"):
                    db_name = os.path.splitext(file_name)[0]
                    db_path = os.path.join(db_dir, db_name)
                    output_file = os.path.join(results_output, f"{db_name}_{query_file_base_name}_results.txt")

                    cmd = f"{blast} -query {query_file_path} -db {db_path} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -out {output_file} -evalue {e_value_threshold}"
                    subprocess.run(cmd, shell=True)

                    with open(output_file, 'r') as f:
                        reader = csv.DictReader(f, fieldnames=fieldnames[:-3], delimiter='\t')
                        for row in reader:
                            if float(row['evalue']) <= e_value_threshold:
                                row['database'] = db_name
                                row['query_file_name'] = query_file_base_name
                                row['file_path'] = output_file
                                results.append(row)

    df = pd.DataFrame(results)

    #Define renamed column headers
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

    if args.report_only_lowest_evalue:
        #Ensure 'evalue' column is numeric
        df['evalue'] = pd.to_numeric(df['evalue'], errors = 'coerce')
        #Use evalue to sort
        df.sort_values(by=['database', 'query_file_name', 'evalue'], ascending=[True, True, True], inplace=True)
        #Reassign df with dropped duplicates
        df = df.drop_duplicates(subset=['database', 'query_file_name'], keep='first')

    df.to_csv(os.path.join(results_output, "results.csv"), index=False)
    print(f"\n \033[92mQueries Complete and stored in {results_output}\033[0m \n")

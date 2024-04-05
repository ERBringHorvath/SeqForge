import os
from Bio import SeqIO
import concurrent.futures

def find_fasta_file(basename, fasta_dir):
    for file in os.listdir(fasta_dir):
        if file.startswith(basename) and file.split('.')[-1] in ['fa', 'fas', 'fasta', 'ffn', 'fna']:
            return os.path.join(fasta_dir, file)
    return None

def process_blast_file(blast_file, blast_dir, fasta_dir, evalue_cutoff):
    file_path = os.path.join(blast_dir, blast_file)
    if os.stat(file_path).st_size == 0:
        print(f"Skipping empty file: {blast_file}")
        return {}

    print(f"\n \033[93mProcessing blast file: {blast_file}\033[0m\n")
    unique_sequences = {}

    base_name = blast_file.split('_')[0]
    original_fasta = find_fasta_file(base_name, fasta_dir)

    if original_fasta is None:
        print(f"\n \033[93mNo matching FASTA file found for {blast_file}. Skipping.\033[0m\n")
        return {}

    with open(os.path.join(blast_dir, blast_file), 'r') as blast_f:
        for line in blast_f:
            parts = line.strip().split('\t')
            if float(parts[10]) < evalue_cutoff:
                seq_id = parts[1]
                sstart = int(parts[8])  # Subject start
                send = int(parts[9])  # Subject end
                # Adjust for the strand; TBLASTN results can span both strands
                strand = 1 if sstart < send else -1

                if seq_id not in unique_sequences:
                    for seq_record in SeqIO.parse(original_fasta, "fasta"):
                        if seq_record.id == seq_id:
                            #Extract based on the strand
                            if strand == 1:
                                aligned_sequence = seq_record.seq[sstart-1:send]
                            else:  #Negative strand, reverse complement
                                aligned_sequence = seq_record.seq[send-1:sstart].reverse_complement()
                            
                            new_record = SeqIO.SeqRecord(aligned_sequence, id=seq_record.id, description="aligned region")
                            unique_sequences[seq_id] = new_record
                            break  #Found the sequence, no need to continue the loop

    return unique_sequences

def extract_sequences_tabular(args):

    blast_dir = args.results_directory
    fasta_dir = args.fasta_directory
    output_fasta = args.output_fasta
    evalue_cutoff = float(args.evalue) if args.evalue else 0.001
    cores = int(args.threads) if args.threads else 1

    blast_files = [file for file in os.listdir(blast_dir) if file.endswith('.txt')]

    print(f"\n \033[92mTotal BLAST files: {len(blast_files)}\033[0m")
    
    all_unique_sequences = {}

    with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
        futures = {executor.submit(process_blast_file, blast_file, blast_dir, fasta_dir, evalue_cutoff): blast_file for blast_file in blast_files}

        for future in concurrent.futures.as_completed(futures):
            unique_sequences = future.result()
            for seq_id, seq_record in unique_sequences.items():
                if seq_id not in all_unique_sequences:
                    all_unique_sequences[seq_id] = seq_record

    print(f"Writing {len(all_unique_sequences)} sequences to output.")      

    SeqIO.write(all_unique_sequences.values(), output_fasta, "fasta")
    print("Sequences processed and written to:", output_fasta)
import os, logging
from Bio import SeqIO
from tqdm import tqdm

def run_split(args):
    input_file = args.input
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    allowed_extensions = ['.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa']

    _, file_extension = os.path.splitext(input_file)
    if file_extension.lower() not in allowed_extensions:
        logging.warning(f"Skipped {input_file} due to incorrect file extension.")
        return
    
    try:
        total_records = sum(1 for _ in SeqIO.parse(input_file, "fasta"))
    except Exception as e:
        logging.error(f"Error parsing {input_file}: {e}")
        return
    
    with tqdm(total=total_records, desc="Processing") as pbar:
        for record in SeqIO.parse(input_file, "fasta"):
            try:
                output_file = os.path.join(output_dir, f"{record.id}.fasta")
                SeqIO.write(record, output_file, "fasta")
                pbar.update(1)
            except Exception as e:
                logging.error(f"Error processing record {record.id} in {input_file}: {e}")
                continue
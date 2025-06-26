import os, logging, gzip
from io import TextIOWrapper
from Bio import SeqIO
from tqdm import tqdm
from math import ceil

def run_split(args):
    input_file = args.input
    output_dir = args.output_dir
    fragment_size = args.fragment
    compress = args.compress

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    allowed_extensions = ['.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa']
    input_basename = os.path.basename(input_file)
    base_name, _ = os.path.splitext(input_basename)

    _, file_extension = os.path.splitext(input_file)
    if file_extension.lower() not in allowed_extensions:
        logging.warning(f"\n \033[91mSkipped {input_file} due to incorrect file extension.\033[0m")
        return

    try:
        records = list(SeqIO.parse(input_file, "fasta"))
    except Exception as e:
        logging.error(f"\n \033[91mError parsing {input_file}: {e}\033[0m")
        return

    total_records = len(records)

    if fragment_size:
        num_chunks = ceil(total_records / fragment_size)
        digits = len(str(num_chunks))
        print(f"\n\033[94mSplitting {total_records} sequences into {num_chunks} file(s) of ~{fragment_size} each\033[0m")

        with tqdm(total=total_records, desc="Fragmenting") as pbar:
            for i in range(num_chunks):
                chunk_records = records[i * fragment_size : (i + 1) * fragment_size]
                chunk_num = str(i + 1).zfill(digits)
                ext = ".fa.gz" if compress else ".fa"
                out_name = f"{base_name}_Frag{fragment_size}_{chunk_num}{ext}"
                out_path = os.path.join(output_dir, out_name)
                try:
                    if compress:
                        with gzip.open(out_path, "wt") as f:
                            SeqIO.write(chunk_records, f, "fasta")
                    else:
                        SeqIO.write(chunk_records, out_path, "fasta")
                    pbar.update(len(chunk_records))
                except Exception as e:
                    logging.error(f"\n \033[91mError writing chunk {chunk_num}: {e}\033[0m")
    else:
        print(f"\n\033[94mSplitting into individual FASTA files...\033[0m")
        with tqdm(total=total_records, desc="Processing") as pbar:
            for record in records:
                out_name = f"{record.id}.fasta.gz" if compress else f"{record.id}.fasta"
                out_path = os.path.join(output_dir, out_name)
                try:
                    if compress:
                        with gzip.open(out_path, "wt") as f:
                            SeqIO.write([record], f, "fasta")
                    else:
                        SeqIO.write(record, out_path, "fasta")
                    pbar.update(1)
                except Exception as e:
                    logging.error(f"\n \033[91mError writing {record.id}: {e}\033[0m")

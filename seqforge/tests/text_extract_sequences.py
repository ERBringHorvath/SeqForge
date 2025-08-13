from pathlib import Path
import pandas as pd
from seqforge.extract_sequences import extract_sequences_from_csv

def test_extract_simple_region(tmp_path, capsys):
    #minimal FASTA and CSV
    fasta = tmp_path / "g.fna"
    fasta.write_text(">contig1\nACGTACGTAC\n")
    csv = tmp_path / "all_results.csv"
    rows = [{
        "qseqid":"q1","sseqid":"contig1","pident":95,"length":4,"mismatch":0,"gapopen":0,
        "qstart":1,"qend":4,"sstart":3,"send":6,"evalue":1e-20,"bitscore":50,"qlen":4,"sframe":1,
        "database":"g","query_file_name":"q"
    }]
    pd.DataFrame(rows).to_csv(csv, index=False)

    out = tmp_path / "out.fna"
    extract_sequences_from_csv(
        csv_path=str(csv),
        fasta_input=str(fasta.parent),
        output_fasta=str(out),
        translate=False,
        evalue=1e-5, min_perc=90.0, min_cov=50.0,
        up=0, down=0, keep_temp_files=False, logger=None, threads=1
    )
    s = out.read_text()
    assert ">contig1_g_q_region" in s
    assert "GTAC" in s
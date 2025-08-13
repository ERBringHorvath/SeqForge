from pathlib import Path
import gzip
from seqforge.utils.file_handler import collect_fasta_files

def test_collect_single_and_gz(tmp_path):
    fa = tmp_path / "a.fa"
    fa.write_text(">a\nACGT\n")
    gz = tmp_path / "b.fasta.gz"
    with gzip.open(gz, "wt") as f:
        f.write(">b\nTTTT\n")

    files, tdir = collect_fasta_files(str(fa))
    assert len(files) == 1 and tdir is None

    files2, tdir2 = collect_fasta_files(str(gz))
    assert len(files2) == 1 and tdir2 is None
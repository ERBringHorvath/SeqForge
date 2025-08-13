from types import SimpleNamespace
from pathlib import Path
import random
from seqforge.generate_unique_fasta_headers import run_unique_fasta_headers

def test_unique_headers_appends_and_random(tmp_path, monkeypatch):
    #deterministic randomness
    monkeypatch.setattr(random, "choices", lambda *a, **k: list("ABCDE"))

    fasta = tmp_path / "x.faa"
    fasta.write_text(">cds1\nAAAA\n>cds2\nTTTT\n")
    outdir = tmp_path / "out"

    args = SimpleNamespace(
        fasta_directory=str(fasta),
        in_place=False,
        output_dir=str(outdir),
        output=None,
        progress="none",
    )
    run_unique_fasta_headers(args)
    out = (outdir / "x.faa").read_text().splitlines()
    assert out[0].startswith(">cds1_x_ABCDE")
    assert out[2].startswith(">cds2_x_ABCDE") is False
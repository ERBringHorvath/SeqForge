from seqforge.utils.progress import ProgressHandler

def test_progress_basic():
    ph = ProgressHandler(total=3, prefix="T", mode="bar")
    ph.update(1)
    ph.update(2)
    ph.finish()
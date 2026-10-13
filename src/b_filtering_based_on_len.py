"""Length filter for FASTQ files.

Reads are copied byte for byte. Standard 4-line FASTQ (as written by ONT
basecallers) is parsed directly, which is several times faster than Biopython;
multi-line FASTQ falls back to Biopython. Gzipped input (.gz) is accepted.
"""
import gzip
import logging
import os

logger = logging.getLogger(__name__)


class _MultiLineFastq(Exception):
    pass


def _open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _filter_four_line(input_fastq, out, min_length):
    kept = 0
    with _open_text(input_fastq) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not header.startswith("@") or not plus.startswith("+") or not qual:
                raise _MultiLineFastq(header[:50])
            if len(seq.rstrip("\r\n")) >= min_length:
                out.write(header)
                out.write(seq)
                out.write(plus)
                out.write(qual if qual.endswith("\n") else qual + "\n")
                kept += 1
    return kept


def _filter_biopython(input_fastq, out, min_length):
    from Bio import SeqIO
    with _open_text(input_fastq) as fh:
        kept = (rec for rec in SeqIO.parse(fh, "fastq") if len(rec.seq) >= min_length)
        return SeqIO.write(kept, out, "fastq")


def filter_fastq_by_length(input_fastq, output_fastq, min_length):
    """Write reads of `input_fastq` with length >= `min_length` to `output_fastq`.

    Returns the number of reads written. The output is written to a temporary
    file and renamed, so a failed run never leaves a partial FASTQ.
    """
    tmp = f"{output_fastq}.tmp{os.getpid()}"
    try:
        try:
            with open(tmp, "w", buffering=1 << 20) as out:
                kept = _filter_four_line(input_fastq, out, min_length)
        except _MultiLineFastq:
            logger.info("Input is not 4-line FASTQ; using Biopython parser")
            with open(tmp, "w") as out:
                kept = _filter_biopython(input_fastq, out, min_length)
        os.replace(tmp, output_fastq)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    logger.info("Kept %d reads >= %d bp", kept, min_length)
    return kept

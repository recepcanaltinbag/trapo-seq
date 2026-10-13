from Bio import SeqIO


def filter_fastq_by_length(input_fastq, output_fastq, min_length):
    """Write reads of `input_fastq` with length >= `min_length` to `output_fastq`.

    Returns the number of reads written.
    """
    with open(output_fastq, "w") as output_handle:
        kept = (rec for rec in SeqIO.parse(input_fastq, "fastq") if len(rec.seq) >= min_length)
        return SeqIO.write(kept, output_handle, "fastq")

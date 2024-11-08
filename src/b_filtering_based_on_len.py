from Bio import SeqIO

def filter_fastq_by_length(input_fastq, output_fastq, min_length):
    with open(output_fastq, "w") as output_handle:
        # Iterate through the FASTQ file and filter sequences by length
        for record in SeqIO.parse(input_fastq, "fastq"):
            if len(record.seq) >= min_length:
                SeqIO.write(record, output_handle, "fastq")

#
# input_fastq = "input_file.fastq"  # Input FASTQ file path
# output_fastq = "filtered_output.fastq"  # Output filtered FASTQ file path
# min_length = 1000  # Set the length threshold
    
# filter_fastq_by_length(input_fastq, output_fastq, min_length)
# print(f"Filtering complete. Sequences >= {min_length} bp saved to {output_fastq}")

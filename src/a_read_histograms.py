#-------------------------------------
#This code used for visualization of read lengths

#INPUT: .fastq file
#OUTPUT: Graphical outputs

#Recep Can Altınbağ, 24 10 2024, v0.0
#-------------------------------------

import matplotlib.pyplot as plt
import math

def calculate_read_lengths(fastq_file):
    """
    Reads a FASTQ file and calculates the length of each sequence (read).
    
    Parameters:
    fastq_file (str): Path to the input FASTQ file.

    Returns:
    list: A list containing the lengths of all reads in the FASTQ file.
    """
    read_lengths = []
    
    with open(fastq_file, 'r') as f:
        # Iterate through the FASTQ file, skipping unnecessary lines
        for i, line in enumerate(f):
            if i % 4 == 1:  # Sequence line (every 4th line, starting from the second)
                read_lengths.append(len(line.strip()))  # Calculate length of the sequence
                
    return read_lengths

def determine_bin_size(read_lengths):
    """
    Determines the bin size for the histogram based on the number of reads.
    
    Parameters:
    read_lengths (list): A list of read lengths.

    Returns:
    int: The bin size for the histogram.
    """
    # Using the square root rule for bin size
    num_reads = len(read_lengths)
    bin_count = math.ceil(math.sqrt(num_reads))
    
    # Calculate the range of the data (max - min)
    data_range = max(read_lengths) - min(read_lengths)
    
    # Determine bin size based on range and number of bins
    bin_size = math.ceil(data_range / bin_count)
    
    return bin_size

def plot_histogram(fastq_file, output):
    """
    Plots a histogram of the read lengths with a dynamically calculated bin size.
    
    Parameters:
    read_lengths (list): A list of read lengths.
    """
    read_lengths = calculate_read_lengths(fastq_file)

    # Determine bin size dynamically
    bin_size = determine_bin_size(read_lengths)
    
    plt.figure(figsize=(10, 6))
    
    # Plot the histogram with calculated bin size
    plt.hist(read_lengths, bins=range(min(read_lengths), max(read_lengths) + bin_size, bin_size), color='skyblue', edgecolor='black')
    
    # Add labels and title for clarity
    plt.title("Read Length Distribution", fontsize=16)
    plt.xlabel("Read Length", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.savefig(f"{output}.pdf")
    print(f'Figure saved as {output}')
    # Show the plot
    #plt.show()

#----------------------------------------------------
# EXAMPLE USAGE -------------------------------------
# INPUTS
# Path to the FASTQ file (replace with your file path)

#fastq_file = "5bac_4000filtered.fastq"

# OUTPUTS
#output = "Read_len_histogram"
#----------------------------------------------------
    
#Plot the histogram
#plot_histogram(fastq_file, output)


#-------------------------------------
#This code used for mapping of reads to the plasmid and genome with using bash script in /scripts/01_data_prep_from_py.sh

#INPUT: .fastq file
#OUTPUT: bam, sam, txt mapping files

#minimap2 is needed!

#Recep Can Altınbağ, 11 11 2024, v0.1
#-------------------------------------





import subprocess
import os



script_path = os.path.abspath("./scripts/01_data_prep_from_py.sh")  # Adjust path as needed

def main_mapping(data_dir, query_fasta, genome_fasta, force):
    print("Script Path:", script_path)
    query_fasta = os.path.abspath(query_fasta)
    genome_fasta = os.path.abspath(genome_fasta)
    args = ["bash", script_path, query_fasta, genome_fasta]
    # If --force is set, append it to the arguments list
    if force:
        args.append("--force")

    with subprocess.Popen(
        args,  # Pass the arguments to the script
        cwd=data_dir,  # Set the working directory to data_dir
        text=True,
        stderr=subprocess.PIPE,
        env={**os.environ, "PYTHONUNBUFFERED": "1"}  # Unbuffered environment
    ) as process:
        # Read the output line by line
        #for line in process.stdout:
        #    print(line, end="")  # Print each line as it’s produced, without extra newlines

        # Read and print errors if any
        for error_line in process.stderr:
            print("Process:", error_line, end="")
        
    # Check if the process ended with an error
    if process.returncode == 0:
        print("\nScript ran successfully!")
    else:
        print("\nScript encountered an error.")



def old_func():
    # Run the script with the arguments in the specified directory
    result = subprocess.run(
        ["bash", script_path, query_fasta, genome_fasta],  # Pass the arguments to the script
        cwd=data_dir,  # Set the working directory to data_dir
        capture_output=True,
        text=True
    )

    # Check the result and print outputs
    if result.returncode == 0:
        print(f"Script ran successfully in the {data_dir} folder!")
        print("Output:", result.stdout)
    else:
        print("Script encountered an error.")
        print("Error:", result.stderr)
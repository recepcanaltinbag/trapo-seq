import argparse
from src.a_read_histograms import plot_histogram
from src.b_filtering_based_on_len import filter_fastq_by_length
from src.a_data_mapping import main_mapping
from src.e_insert_finder_from_bam import main_insert_finder_from_bam
from src.g_blast_annot import main_annot
from src.h_is_stats import main_is_stat
from src.j_dr_finder import main_dr_finder

trapo_seq_logo_v2 = """
+-+-+-+-+-+-+-+-+-+
|T|R|A|P|O|-|S|E|Q|
+-+-+-+-+-+-+-+-+-+
"""

trapo_seq_logo = """
 ___________________________________________
| ___ ____ ____ ___  ____    ____ ____ ____ |
|  |  |__/ |__| |__] |  | __ [__  |___ |  | |
|  |  |  \ |  | |    |__|    ___] |___ |_\| | 
|___________________________________________|

________Traposome Sequencing Pipeline_______

"""

def main():
    print("\033[93m" + trapo_seq_logo + "\033[0m")

    # Main Parser
    parser = argparse.ArgumentParser(description="Traposome Sequencing Pipeline")
    parser.add_argument('--all-help', action='store_true', help="Show complete help message")

    subparsers = parser.add_subparsers(dest="command")

    # Subcommands
    read_histogram_parser = subparsers.add_parser("read_histogram", help="creates histogram of reads based on lengths")
    read_histogram_parser.add_argument("-v", "--verbose", action="store_true", help="Extented output")
    read_histogram_parser.add_argument("-f", "--fastq", type=str, required=True, help="Path of fastq file")
    read_histogram_parser.add_argument("-o", "--output", type=str, required=True, help="path of output")

    filter_parser = subparsers.add_parser("filter", help="filter reads based on lengths")
    filter_parser.add_argument("-f", "--fastq", type=str, required=True, help="Path of fastq file")
    filter_parser.add_argument("-o", "--output", type=str, required=True, help="path of output filtered fastq file")
    filter_parser.add_argument("-l", "--length", type=int, default=1000, help="len of filter")

    map_parser = subparsers.add_parser("map", help="map reads to plasmid and genome")
    map_parser.add_argument("-d", "--input_dir", type=str, required=True, help="Path of fastq file including folders")
    map_parser.add_argument("-p", "--plasmid", type=str, required=True, help="Path of trap plasmid in fasta format")
    map_parser.add_argument("-g", "--genome", type=str, required=True, help="path of genome in fasta format")
    map_parser.add_argument("-f", "--force", action="store_true", help="Force overwrite if the file exists")

    insert_finder_parser = subparsers.add_parser("insert_finder", help="finding inserts from .bam file")
    insert_finder_parser.add_argument("-b", "--bam", type=str, required=True, help="Path of bam file")
    insert_finder_parser.add_argument("-o", "--output", type=str, required=True, help="Out File")
    insert_finder_parser.add_argument("-t", "--threshold", type=int, default=500, help="Insertion Threshold, insertions lower than this will be disregarded!")

    blast_annot_parser = subparsers.add_parser("blast_annot", help="Annotation with Blast")
    blast_annot_parser.add_argument("-i", "--ins_bam", type=str, required=True, help="Path of insertions from bam tabular file")
    blast_annot_parser.add_argument("-m", "--mapped_fasta", type=str, required=True, help="path of mapped reads in fasta format")
    blast_annot_parser.add_argument("-g", "--genome_fasta", type=str, required=True, help="path of genome in fasta format")
    blast_annot_parser.add_argument("--is_fasta", type=str, required=True, help="path of manually curated ISes")
    blast_annot_parser.add_argument("--temp", type=str, default="data/temp", help="temp directory, default :temp")
    blast_annot_parser.add_argument("-o", "--output", type=str, default="best_alignments.tab", help="Output file")
    blast_annot_parser.add_argument("-t", "--threshold", type=int, default=70, help="Threshold for elimination of Tns, If lower than this, disregard. Value: 0-100, default:70")
    blast_annot_parser.add_argument("--threads", type=int, default=2, help="How many threads it is wanted for blast searches, default 2")
    blast_annot_parser.add_argument("--debug", action="store_true", help="For debugging, it can be given. Prints everything")
    blast_annot_parser.add_argument("--no_temp", action="store_false", help="To delete temp folder, add this")
    blast_annot_parser.add_argument("--partial_threshold", type=int, default=80, help="For big insertions, if it has lower coverage than this value (and partial_len is used with this value), eliminate. Value: 0-100, default:80")
    blast_annot_parser.add_argument("--partial_len", type=int, default=8000, help="For big insertions, if it is higher than this value, elimiation will be based on partial threshold, It is effective if there is multiple IS and TNs together. default:8000")

    is_stat_parser = subparsers.add_parser("is_stat", help="Stats of ISes and Tn's")
    is_stat_parser.add_argument("-d", "--input_dir", type=str, required=True, help="Path of data folder")
    is_stat_parser.add_argument("-o", "--output", type=str, required=True, help="path of output stat .rcp file")

    heatmap_parser = subparsers.add_parser("heatmap", help="Creates a heatmap from *.rcp stats file")
    heatmap_parser.add_argument("-r", "--rcp_file", type=str, required=True, help="Path of *.rcp stats file")
    heatmap_parser.add_argument("-o", "--output", type=str, required=True, help="Output file")
    heatmap_parser.add_argument("-t", "--out_text", type=str, required=True, help="Output file")

    dr_finder_parser = subparsers.add_parser("dr_finder", help="Finding DRs")
    dr_finder_parser.add_argument("-d", "--input_dir", type=str, required=True, help="Path of data folder")
    dr_finder_parser.add_argument("-p", "--plasmid", type=str, required=True, help="Path of trap plasmid in fasta format")
    dr_finder_parser.add_argument("-g", "--gap", type=int, default=500, help="Gap, to extract more left and right flanking sequences from reads, deafult: 500")


    # Arguments
    args = parser.parse_args()

    #A_READ_HISTOGRAMS
    if args.all_help:
        parser.print_help()
        print("\n" + "-"*40 + "\n")
        for subcommand_name, subcommand_parser in subparsers.choices.items():
            print(f"\n\033[93mSubcommand: {subcommand_name}\033[0m")
            subcommand_parser.print_help()
            print("\n" + "-"*40 + "\n")
    else:
        if args.command == "read_histogram":
            if args.verbose:
                print('Histograms')
            else:
                print("Histograms...")
                plot_histogram(args.fastq, args.output)
        #B_FILTER
        elif args.command == "filter":
            print(f"Filter process, Input: {args.fastq}, Output: {args.output}, Filter Len: {args.length}")
            filter_fastq_by_length(args.fastq, args.output, args.length)
            print('Filter process was ended')          
        #
        #C_ALIGNMENT
        elif args.command == "map":
            print(f"\nMapping process\nInput Directory: {args.input_dir}, \nPlasmid: {args.plasmid}, \nGenome: {args.genome}, \nOverwriting: {args.force}")
            main_mapping(args.input_dir, args.plasmid, args.genome, args.force)
            print('Mapping process was ended')

        #E_INSERT_FINDER_FROM_BAM
        elif args.command == "insert_finder":
            print(f"\nMapping process\nInput BAM: {args.bam}, \nOutput: {args.output}, \nThreshold: {args.threshold}")
            main_insert_finder_from_bam(args.bam, args.output, args.threshold)


        #G_BLAST_THE_INSERTS_
        elif args.command == "blast_annot":
            print(f"\nBlasting to annotate insertions..\n")
            main_annot(args.ins_bam, args.mapped_fasta, args.genome_fasta, args.is_fasta, args.temp, args.output, args.threshold, args.threads, args.debug, args.no_temp, args.partial_threshold, args.partial_len)

        #H_IS_Stats
        elif args.command == "is_stat":
            print(f"\nStats of ISes..\n")
            main_is_stat(args.input_dir, args.output)

        #I_Heatmap
        elif args.command == "heatmap":
            from src.i_heatmap import main_heatmap
            print(f"\nCreating Heatmap..\n")
            main_heatmap(args.rcp_file, args.output, args.out_text)

        #J_DR_finder
        elif args.command == "dr_finder":
            print(f"\nFinding DRs..\n")
            main_dr_finder(args.input_dir, args.plasmid, args.gap)

        else:
            print("Please enter a command.")


if __name__ == "__main__":
    main()
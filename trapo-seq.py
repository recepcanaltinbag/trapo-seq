import argparse
import logging
import subprocess
import os
import sys

REPO_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO_DIR)

from src import __version__

trapo_seq_logo_v2 = """
+-+-+-+-+-+-+-+-+-+
|T|R|A|P|O|-|S|E|Q|
+-+-+-+-+-+-+-+-+-+
"""

trapo_seq_logo = r"""
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
    parser.add_argument(
        '--version',
        action='version',
        version=f'trapo-seq {__version__}',
        help="Show program's version number and exit"
    )


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
    map_parser.add_argument("-t", "--threads", type=int, default=4, help="Threads for minimap2 and samtools, default: 4")
    map_parser.add_argument("--linear", action="store_true", help="Map to the plasmid as a linear sequence only (old behaviour); by default the circular plasmid is handled")

    insert_finder_parser = subparsers.add_parser("insert_finder", help="finding inserts from .bam file")
    insert_finder_parser.add_argument("-b", "--bam", type=str, required=True, help="Path of bam file")
    insert_finder_parser.add_argument("-o", "--output", type=str, required=True, help="Out File")
    insert_finder_parser.add_argument("-t", "--threshold", type=int, default=500, help="Insertion Threshold, insertions lower than this will be disregarded!")

    insert_finder_batch_parser = subparsers.add_parser("insert_finder_batch", help="finding inserts from .bam file in Batch mode")
    insert_finder_batch_parser.add_argument("-d", "--input_dir", type=str, required=True, help="Folder of folders including bam files")
    insert_finder_batch_parser.add_argument("-t", "--threshold", type=int, default=500, help="Insertion Threshold, insertions lower than this will be disregarded!")

    blast_annot_parser = subparsers.add_parser("blast_annot", help="Annotation with Blast")
    blast_annot_parser.add_argument("-i", "--ins_bam", type=str, required=True, help="Path of insertions from bam tabular file")
    blast_annot_parser.add_argument("-m", "--mapped_fasta", type=str, required=True, help="path of mapped reads in fasta format")
    blast_annot_parser.add_argument("-g", "--genome_fasta", type=str, required=True, help="path of genome in fasta format")
    blast_annot_parser.add_argument("--is_fasta", type=str, default=None, help="path of manually curated ISes (>Family_Name headers); a placeholder database is used if omitted")
    blast_annot_parser.add_argument("--temp", type=str, default="temp", help="temp directory, default: temp")
    blast_annot_parser.add_argument("-o", "--output", type=str, default="best_alignment.tab", help="Output file")
    blast_annot_parser.add_argument("-t", "--threshold", type=int, default=70, help="Threshold for elimination of Tns, If lower than this, disregard. Value: 0-100, default:70")
    blast_annot_parser.add_argument("--threads", type=int, default=2, help="How many threads it is wanted for blast searches, default 2")
    blast_annot_parser.add_argument("--debug", action="store_true", help="For debugging, it can be given. Prints everything")
    blast_annot_parser.add_argument("--keep_temp", action="store_true", help="Keep temporary BLAST files and databases (for debugging)")
    blast_annot_parser.add_argument("--evalue", type=float, default=1e-10, help="BLAST E-value cutoff for insert annotation, default: 1e-10")
    blast_annot_parser.add_argument("--blast_task", choices=["megablast", "dc-megablast", "blastn"], default="megablast", help="BLAST search mode; dc-megablast is more sensitive for divergent inserts, default: megablast")
    blast_annot_parser.add_argument("--partial_threshold", type=int, default=80, help="For big insertions, if it has lower coverage than this value (and partial_len is used with this value), eliminate. Value: 0-100, default:80")
    blast_annot_parser.add_argument("--partial_len", type=int, default=8000, help="For big insertions, if it is higher than this value, elimiation will be based on partial threshold, It is effective if there is multiple IS and TNs together. default:8000")

    blast_annot_batch_parser = subparsers.add_parser("blast_annot_batch", help="Annotation with Blast in Batch mode")
    blast_annot_batch_parser.add_argument("-d", "--input_dir", type=str, required=True, help="Folder of folders including bam files")
    blast_annot_batch_parser.add_argument("-g", "--genome_fasta", type=str, required=True, help="path of genome in fasta format")
    blast_annot_batch_parser.add_argument("--is_fasta", type=str, default=None, help="path of manually curated ISes (>Family_Name headers); a placeholder database is used if omitted")
    blast_annot_batch_parser.add_argument("--temp", type=str, default=None, help="temp directory, default: <input_dir>/temp")
    blast_annot_batch_parser.add_argument("-t", "--threshold", type=int, default=70, help="Threshold for elimination of Tns, If lower than this, disregard. Value: 0-100, default:70")
    blast_annot_batch_parser.add_argument("--threads", type=int, default=2, help="How many threads it is wanted for blast searches, default 2")
    blast_annot_batch_parser.add_argument("--debug", action="store_true", help="For debugging, it can be given. Prints everything")
    blast_annot_batch_parser.add_argument("--keep_temp", action="store_true", help="Keep temporary BLAST files and databases (for debugging)")
    blast_annot_batch_parser.add_argument("--evalue", type=float, default=1e-10, help="BLAST E-value cutoff for insert annotation, default: 1e-10")
    blast_annot_batch_parser.add_argument("--blast_task", choices=["megablast", "dc-megablast", "blastn"], default="megablast", help="BLAST search mode; dc-megablast is more sensitive for divergent inserts, default: megablast")
    blast_annot_batch_parser.add_argument("--partial_threshold", type=int, default=80, help="For big insertions, if it has lower coverage than this value (and partial_len is used with this value), eliminate. Value: 0-100, default:80")
    blast_annot_batch_parser.add_argument("--partial_len", type=int, default=8000, help="For big insertions, if it is higher than this value, elimiation will be based on partial threshold, It is effective if there is multiple IS and TNs together. default:8000")

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
    dr_finder_parser.add_argument("-g", "--gap", type=int, default=500, help="Flank bases taken from the read on each side of the insertion, default: 500")
    dr_finder_parser.add_argument("-t", "--threshold", type=int, default=200, help="Max distance (bp, inclusive) between a plasmid flank and the insertion edge on the read, default: 200")
    dr_finder_parser.add_argument("-r", "--repeat_threshold", type=int, default=30, help="Max DR length (inclusive), default: 30")
    dr_finder_parser.add_argument("--is_fasta", type=str, default=None, help="Curated IS FASTA; IS lengths sharpen DR boundaries for IS_DB insertions")
    dr_finder_parser.add_argument("--threads", type=int, default=2, help="BLAST threads, default: 2")

    dr_logo_parser = subparsers.add_parser("dr_logo", help="Finding DR Logos")
    dr_logo_parser.add_argument("-d", "--input_dir", type=str, required=True, help="Path of data folder")
    dr_logo_parser.add_argument("-p", "--plasmid", type=str, required=True, help="Path of trap plasmid in fasta format")
    dr_logo_parser.add_argument("-o", "--out_dir", type=str, required=True, help="Output folder")
    dr_logo_parser.add_argument("-g", "--gap_threshold", type=int, default=60, help="Max %% gaps per alignment column kept in the logo (inclusive), default: 60")

    in_del_plot_parser = subparsers.add_parser("in_del_plot", help="Plotting indels and dels to analyze possible excisions and similar trends")
    in_del_plot_parser.add_argument("-i", "--ins_bam", type=str, required=True, help="Path of insertions from bam tabular file")
    in_del_plot_parser.add_argument("-b", "--bam", type=str, required=True, help="Path of bam file")
    in_del_plot_parser.add_argument("-o", "--output", type=str, required=True, help="Out File")
    in_del_plot_parser.add_argument("--in_threshold", type=int, default=2, help="Disregard lower insertion than this, defult: 2")
    in_del_plot_parser.add_argument("--del_threshold", type=int, default=5, help="Disregard lower deletions than this, defult: 5")

    kde_mobile_parser = subparsers.add_parser("kde_mobile", help="Plotting kdes, to anaylze jumping DNA sizes")
    kde_mobile_parser.add_argument("-b", "--bam", type=str, required=True, help="Path of bam file")
    kde_mobile_parser.add_argument("-o", "--output", type=str, required=True, help="Out File")
    kde_mobile_parser.add_argument("--plasmid_threshold", type=int, default=2000, help="Plasmid len threshold, if the alignment len is lower than this disregard")

    map_dist_parser = subparsers.add_parser("map_dist", help="Plotting dist, to anaylze jumping DNA coordinates both genome and plasmid")
    map_dist_parser.add_argument("-b", "--bam", type=str, required=True, help="Path of bam file")
    map_dist_parser.add_argument("-o", "--output", type=str, required=True, help="Out File")
    map_dist_parser.add_argument("-f", "--fastq", type=str, default=None, help="Optional raw fastq, to report the fraction of reads mapped")

    novel_parser = subparsers.add_parser("novel_elements", help="Cluster and characterise novel inserts (Genome and Unannotated): boundaries, TIRs, DRs, longest ORF, read consensus. Run dr_finder first for plasmid DRs; needs minimap2 and mafft")
    novel_parser.add_argument("-d", "--input_dir", type=str, required=True, help="Data folder with barcode folders (after blast_annot_batch)")
    novel_parser.add_argument("-g", "--genome_fasta", type=str, required=True, help="Genome FASTA, the same one used for blast_annot_batch")
    novel_parser.add_argument("--is_fasta", type=str, default=None, help="Curated IS FASTA, to flag clusters that are copies of known IS (recommended)")
    novel_parser.add_argument("-o", "--out_dir", type=str, default=None, help="Output folder, default: <input_dir>/novel_elements")
    novel_parser.add_argument("--min_reads", type=int, default=2, help="Minimum reads per element cluster, default: 2")
    novel_parser.add_argument("--threads", type=int, default=2, help="Threads for minimap2, default: 2")
    novel_parser.add_argument("--flank_search", type=int, default=60, help="bp searched around each estimated element end for TIRs, default: 60")
    novel_parser.add_argument("--min_tir_len", type=int, default=10, help="Minimum TIR length, default: 10")
    novel_parser.add_argument("--min_tir_identity", type=float, default=0.75, help="Minimum TIR identity (0-1), default: 0.75")
    novel_parser.add_argument("--min_tir_score", type=int, default=22, help="Minimum TIR alignment score, default: 22")
    novel_parser.add_argument("--min_dr_len", type=int, default=3, help="Minimum genomic DR length counted as evidence, default: 3")
    novel_parser.add_argument("--min_orf_aa", type=int, default=200, help="Minimum longest ORF (aa) counted as evidence, default: 200")
    novel_parser.add_argument("--keep_temp", action="store_true", help="Keep temporary files")
    novel_parser.add_argument("--no_unannotated", action="store_true", help="Skip read-consensus clusters of Unannotated inserts")
    novel_parser.add_argument("--max_consensus_reads", type=int, default=30, help="Reads used per consensus, default: 30")
    novel_parser.add_argument("--consensus_flank", type=int, default=20, help="bp searched for TIRs at consensus ends, default: 20")
    novel_parser.add_argument("--consensus_base_fraction", type=float, default=0.4, help="Min fraction of reads supporting a consensus base, default: 0.4")
    novel_parser.add_argument("--polish_rounds", type=int, default=2, help="Consensus polishing rounds, default: 2")
    novel_parser.add_argument("--dr_dir", type=str, default=None, help="dr_finder output folder, default: <input_dir>/insertions")

    test_parser = subparsers.add_parser("test", help="to test the pipeline")

    # Arguments
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

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
            from src.a_read_histograms import plot_histogram
            print("Histograms...")
            plot_histogram(args.fastq, args.output)
        
        #B_FILTER
        elif args.command == "filter":
            from src.b_filtering_based_on_len import filter_fastq_by_length
            print(f"Filter process, Input: {args.fastq}, Output: {args.output}, Filter Len: {args.length}")
            filter_fastq_by_length(args.fastq, args.output, args.length)
            print('Filter process was ended')          
        
        #C_ALIGNMENT
        elif args.command == "map":
            from src.a_data_mapping import main_mapping
            print(f"\nMapping process\nInput Directory: {args.input_dir}, \nPlasmid: {args.plasmid}, \nGenome: {args.genome}, \nOverwriting: {args.force}")
            from src.a_data_mapping import MappingError
            try:
                main_mapping(args.input_dir, args.plasmid, args.genome, args.force, threads=args.threads,
                             circular=not args.linear)
            except (MappingError, FileNotFoundError) as err:
                sys.exit(f"Mapping failed: {err}")
            print('Mapping process was ended')

        #E_INSERT_FINDER_FROM_BAM
        elif args.command == "insert_finder":
            from src.e_insert_finder_from_bam import main_insert_finder_from_bam
            print(f"\nMapping process\nInput BAM: {args.bam}, \nOutput: {args.output}, \nThreshold: {args.threshold}")
            main_insert_finder_from_bam(args.bam, args.output, args.threshold)

        #BATCH VERSION
        elif args.command == "insert_finder_batch":
            from src.e_insert_finder_from_bam_batch import main_insert_finder_from_bam_batch
            print(f"\n[BATCH] Mapping process\nInput BAM: {args.input_dir}, \nThreshold: {args.threshold}")
            main_insert_finder_from_bam_batch(args.input_dir, args.threshold)

        #G_BLAST_THE_INSERTS_
        elif args.command == "blast_annot":
            from src.g_blast_annot import main_annot
            print(f"\nBlasting to annotate insertions..\n")
            main_annot(args.ins_bam, args.mapped_fasta, args.genome_fasta, is_fasta=args.is_fasta, temp_dir=args.temp,
                       output_csv=args.output, threshold=args.threshold, threads=args.threads, debug=args.debug,
                       keep_temp=args.keep_temp, partial_threshold=args.partial_threshold, partial_len=args.partial_len,
                       evalue=args.evalue, task=args.blast_task)
        
        #G_BLAST_ANNOT_BATCH
        elif args.command == "blast_annot_batch":
            from src.g_blast_annot_batch import main_annot_batch
            print(f"\n[BATCH] Blasting to annotate insertions..\n")
            main_annot_batch(args.input_dir, args.genome_fasta, is_fasta=args.is_fasta, temp_dir=args.temp,
                             threshold=args.threshold, threads=args.threads, debug=args.debug,
                             keep_temp=args.keep_temp, partial_threshold=args.partial_threshold,
                             partial_len=args.partial_len, evalue=args.evalue, task=args.blast_task)

        #H_IS_Stats
        elif args.command == "is_stat":
            from src.h_is_stats import main_is_stat
            print(f"\nStats of ISes..\n")
            main_is_stat(args.input_dir, args.output)

        #I_Heatmap
        elif args.command == "heatmap":
            from src.i_heatmap import main_heatmap
            print(f"\nCreating Heatmap..\n")
            main_heatmap(args.rcp_file, args.output, args.out_text)

        #J_DR_finder
        elif args.command == "dr_finder":
            from src.j_dr_finder import main_dr_finder
            print(f"\nFinding DRs..\n")
            main_dr_finder(args.input_dir, args.plasmid, args.gap, args.threshold, args.repeat_threshold,
                           threads=args.threads, is_fasta=args.is_fasta)

        #K_DR_Logo
        elif args.command == "dr_logo":
            from src.k_dr_logo import main_dr_logo
            print(f"\nFinding DR Logos..\n")
            main_dr_logo(args.plasmid, args.input_dir, args.out_dir, args.gap_threshold)

        #F_INDEL_PLOTS
        elif args.command == "in_del_plot":
            from src.f_in_del_plot import main_in_del_plot
            print("\nPlotting indels...\n")
            main_in_del_plot(args.ins_bam, args.bam, args.output, args.in_threshold, args.del_threshold)
        
        #D_KDE
        elif args.command == "kde_mobile":
            from src.d_kde_mobile import main_kde_mobile
            print("\nPlotting kdes of mobile bpairs of genome...\n")
            main_kde_mobile(args.bam, args.output, args.plasmid_threshold)
        
        #C_MAP_DIST
        elif args.command == "map_dist":
            from src.c_map_dist import main_map_dist
            print("\nPlotting plasmid and genome read distributions...\n")
            main_map_dist(args.bam, args.output, fastq_file=args.fastq)

        #L_NOVEL_ELEMENTS
        elif args.command == "novel_elements":
            from src.l_novel_elements import main_novel_elements
            print("\nCharacterising novel elements...\n")
            main_novel_elements(args.input_dir, args.genome_fasta, out_dir=args.out_dir,
                                min_reads=args.min_reads, threads=args.threads,
                                flank_search=args.flank_search, is_fasta=args.is_fasta,
                                keep_temp=args.keep_temp, min_tir_len=args.min_tir_len,
                                min_tir_identity=args.min_tir_identity, min_tir_score=args.min_tir_score,
                                min_dr_len=args.min_dr_len, min_orf_aa=args.min_orf_aa,
                                unannotated=not args.no_unannotated, max_consensus_reads=args.max_consensus_reads,
                                consensus_flank=args.consensus_flank, dr_dir=args.dr_dir,
                                consensus_base_fraction=args.consensus_base_fraction, polish_rounds=args.polish_rounds)

        #TEST
        elif args.command == "test":
            print('Starting test!\n')
            subprocess.run(["bash", os.path.join(REPO_DIR, "test_script.sh")], check=True)


        else:
            print("Please enter a command.")


if __name__ == "__main__":
    main()
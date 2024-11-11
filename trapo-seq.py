import argparse
from src.a_read_histograms import plot_histogram
from src.b_filtering_based_on_len import filter_fastq_by_length
from src.a_data_mapping import main_mapping

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

    # Ana parser oluşturma
    parser = argparse.ArgumentParser(description="Traposome Sequencing Pipeline")
    subparsers = parser.add_subparsers(dest="command")

    # islem1 alt komutunu ekleme
    read_histogram_parser = subparsers.add_parser("read_histogram", help="creates histogram of reads based on lengths")
    read_histogram_parser.add_argument("-v", "--verbose", action="store_true", help="Extented output")
    read_histogram_parser.add_argument("-f", "--fastq", type=str, required=True, help="Path of fastq file")
    read_histogram_parser.add_argument("-o", "--output", type=str, required=True, help="path of output")


    filter_parser = subparsers.add_parser("filter", help="filter reads based on lengths")
    filter_parser.add_argument("-f", "--fastq", type=str, required=True, help="Path of fastq file")
    filter_parser.add_argument("-o", "--output", type=str, required=True, help="path of output filtered fastq file")
    filter_parser.add_argument("-l", "--length", type=int, default=1000, help="len of filter")


    filter_parser = subparsers.add_parser("map", help="map reads to plasmid and genome")
    filter_parser.add_argument("-d", "--input_dir", type=str, required=True, help="Path of fastq file including folders")
    filter_parser.add_argument("-p", "--plasmid", type=str, required=True, help="Path of trap plasmid in fasta format")
    filter_parser.add_argument("-g", "--genome", type=str, required=True, help="path of genome in fasta format")
    filter_parser.add_argument("-f", "--force", action="store_true", help="Force overwrite if the file exists")







    # islem2 alt komutunu ekleme
    islem2_parser = subparsers.add_parser("islem2", help="İşlem 2'yi çalıştırır")
    islem2_parser.add_argument("-p", "--parametre", type=str, required=True, help="İşlem 2 için gerekli parametre")
    islem2_parser.add_argument("-n", "--number", type=int, default=1, help="İşlem tekrarı için sayı belirler (varsayılan: 1)")

    # Argümanları ayrıştırma
    args = parser.parse_args()

    #A_READ_HISTOGRAMS
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
        print(f"\nMapping process\n, Input Directory: {args.input_dir}, \nPlasmid: {args.plasmid}, \nGenome: {args.genome}, \nOverwriting: {args.force}")
        main_mapping(args.input_dir, args.plasmid, args.genome, args.force)
        print('Mapping process was ended')
























    elif args.command == "islem2":
        print(f"İşlem 2 çalıştırılıyor, parametre: {args.parametre}, tekrar sayısı: {args.number}")
    else:
        print("Lütfen bir komut belirtin.")


if __name__ == "__main__":
    main()
import argparse
from src.a_read_histograms import plot_histogram


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
    islem1_parser = subparsers.add_parser("read_histogram", help="creates histogram of reads based on lengths")
    islem1_parser.add_argument("-v", "--verbose", action="store_true", help="Extented output")
    islem1_parser.add_argument("-f", "--fastq", type=str, required=True, help="Path of fastq file")
    islem1_parser.add_argument("-o", "--output", type=str, required=True, help="path of output")


    # islem2 alt komutunu ekleme
    islem2_parser = subparsers.add_parser("islem2", help="İşlem 2'yi çalıştırır")
    islem2_parser.add_argument("-p", "--parametre", type=str, required=True, help="İşlem 2 için gerekli parametre")
    islem2_parser.add_argument("-n", "--number", type=int, default=1, help="İşlem tekrarı için sayı belirler (varsayılan: 1)")

    # Argümanları ayrıştırma
    args = parser.parse_args()

    # Komutları çalıştırma
    if args.command == "read_histogram":
        if args.verbose:
            print('Histograms')
        else:
            print("Histograms...")
            plot_histogram(args.fastq, args.output)
    elif args.command == "islem2":
        print(f"İşlem 2 çalıştırılıyor, parametre: {args.parametre}, tekrar sayısı: {args.number}")
    else:
        print("Lütfen bir komut belirtin.")


if __name__ == "__main__":
    main()
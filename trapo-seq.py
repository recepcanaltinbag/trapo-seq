import argparse

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
    islem1_parser = subparsers.add_parser("islem1", help="İşlem 1'i çalıştırır")
    islem1_parser.add_argument("-v", "--verbose", action="store_true", help="Detaylı çıktı sağlar")

    # islem2 alt komutunu ekleme
    islem2_parser = subparsers.add_parser("islem2", help="İşlem 2'yi çalıştırır")
    islem2_parser.add_argument("-p", "--parametre", type=str, required=True, help="İşlem 2 için gerekli parametre")
    islem2_parser.add_argument("-n", "--number", type=int, default=1, help="İşlem tekrarı için sayı belirler (varsayılan: 1)")

    # Argümanları ayrıştırma
    args = parser.parse_args()

    # Komutları çalıştırma
    if args.command == "islem1":
        if args.verbose:
            print("İşlem 1 detaylı modda çalıştırılıyor...")
        else:
            print("İşlem 1 çalıştırılıyor...")
    elif args.command == "islem2":
        print(f"İşlem 2 çalıştırılıyor, parametre: {args.parametre}, tekrar sayısı: {args.number}")
    else:
        print("Lütfen bir komut belirtin.")


if __name__ == "__main__":
    main()
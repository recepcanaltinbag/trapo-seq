import numpy as np
import matplotlib.pyplot as plt

# Genom verisi
'''
genome_length = 6_000_000  # 6 Mbp'lik genom
locations = np.linspace(0, genome_length, 1000)  # 1000 eşit bölge
counts = np.random.randint(0, 200, size=len(locations))  # Rastgele count değerleri (daha düşük yoğunluk)
print(locations)
'''


# Genom verisi
'''
genome_length = 6_000_000  # 6 Mbp'lik genom
locations = np.linspace(0, genome_length, 1000)  # 1000 eşit bölge
counts = np.random.randint(0, 200, size=len(locations))  # Rastgele count değerleri (daha düşük yoğunluk)
print(locations)
'''

def main_polar(genome_array):
    genome_length = len(genome_array)
    counts = []
    locations = np.linspace(0, genome_length, 100)  # 500 eşit bölge
    # Bölgelere göre işlem yapma
    for i in range(1, len(locations)):
        start_idx = int(locations[i-1])
        end_idx = int(locations[i])
        
        # Burada dilimleme yapılıyor: genome_array[start_idx:end_idx]
        counts.append(max(genome_array[start_idx:end_idx]))

    counts.append(0)
    counts = np.array(counts)
    # Dairesel koordinatlar
    theta = 2 * np.pi * locations / genome_length  # Açı hesaplama

    # Renk skalası (yoğunluk için)
    colors = plt.cm.viridis(counts / max(counts))  # Viridis renk haritası

    # En yoğun 5 nokta (count değerleri en yüksek olanlar)
    top_5_indices = np.argsort(counts)[-5:]  # En yüksek 5 değeri tespit et
    top_5_locations = locations[top_5_indices][::-1]  # En yoğun 5 noktanın koordinatları
    top_5_counts = counts[top_5_indices][::-1]  # En yoğun 5 noktanın count değerleri
    top_5_theta = theta[top_5_indices][::-1]  # En yoğun 5 noktanın theta (açı) değerleri

    # Polar plot
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(12, 6))
    # Polar eksende 0 noktasını en tepeye yerleştirme
    ax.set_theta_offset(np.pi / 2)  # 0'ı tepeye koymak için (90 derece)
    ax.set_theta_direction(-1)  # Saat yönünde artan açı
    # Her bir peak için çizim (merkeze doğru yoğunluk artacak şekilde)
    for i in range(len(locations)):
        ax.plot([theta[i], theta[i]], [max(counts) - counts[i], max(counts)], color=colors[i], linewidth=1.5)

    # Genom koordinatlarını dairenin çevresine yerleştirme
    # Bunu her 1000 base pair'lik bölgeye yerleştiriyoruz, yani her 1000 bp'de bir etiket olacak
    num_ticks = 10  # 6 etiket olacak
    tick_positions = np.linspace(0, genome_length, num_ticks)
    tick_labels = [f"{int(pos / 1e4)} Mb" for pos in tick_positions]


    # Tam 1 MB'lik aralıklarla tickler
    num_ticks = 6  # 6 etiket olacak (0, 1, 2, 3, 4, 5 MB)
    tick_positions = np.arange(0, genome_length, 1_000_000)  # 1 MB'lik adımlarla
    tick_labels = [f"{int(pos / 1e6)} Mb" for pos in tick_positions]  # Etiketler (MB cinsinden)


    # Polar grafikteki açılar ile genom bölgesindeki konumları eşleştirmek
    theta_ticks = 2 * np.pi * tick_positions / genome_length

    ax.set_xticks(theta_ticks)  # Dairedeki doğru açılara etiketleri yerleştir
    ax.set_xticklabels(tick_labels)  # Etiketleri yerleştir

    for i in theta_ticks:
        ax.plot([i, i], [0, max(counts)], color='black', lw=1) 

    # En yoğun 5 nokta için işaretler
    for i in range(len(top_5_locations)):
        ax.plot(top_5_theta[i], max(counts), 'bo', markersize=10)  # Yoğun noktaları mavi ile işaretle
        
        # Etiketleri daire dışına yatay olarak yerleştir
        ax.text(top_5_theta[i], max(counts) - max(counts)/10, f"P{i+1}",
                horizontalalignment='center', verticalalignment='center', 
                fontsize=10, color='black', rotation=top_5_theta[i],
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))  # Yatay olarak etiket yerleştir

    # Etiketler için yeni bir alan ekleyelim
    ax2 = fig.add_axes([0.85, 0.1, 0.1, 0.8])  # Bu alan grafiğin yan tarafında olacak
    ax2.axis('off')  # Grafik çizilmeyecek, sadece etiketler olacak

    # Etiketler
    #ax2.text(0, 0.95, "En Yoğun 5 Peak", fontsize=14, fontweight='bold')
    for i in range(len(top_5_locations)):
        ax2.text(0, 0.9 - i * 0.07, f"{i+1}. Peak: {int(top_5_locations[i] / 1e3)} kbp, Count: {top_5_counts[i]}", 
                fontsize=12, ha='center')

    # Görsel detaylar
    ax.set_title("Genom Koordinatlarına Göre Dairesel Peak Grafiği (Yatay Etiketler)", va='bottom')
    ax.set_ylim(0, max(counts))  # Y ekseni sınırları

    ax.yaxis.set_visible(False)  # Merkezden dışa doğru yoğunluk kılavuz çizgileri (azalan sıralama)
    # En yüksek değeri merkeze koy
    guide_values = np.linspace(0, max(counts), 5)  # Kılavuz çizgileri için değerler (0'dan en yüksek değere kadar)
    # Diğer kılavuz çizgileri
    for gv in guide_values[:-1]:  # En yüksek değeri tekrar etmemek için '[:-1]' kısmı
        ax.plot(np.linspace(0, 2 * np.pi, 500), [gv] * 500, '--', color='gray', alpha=0.5, linewidth=0.5)

    sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(0, max(counts)))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='horizontal', pad=0.1, shrink=0.8)
    cbar.set_label("Count Yoğunluğu")

    plt.show()

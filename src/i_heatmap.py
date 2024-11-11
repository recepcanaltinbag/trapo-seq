import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
import os

# Pastel teal renk paletini oluşturma
def create_pastel_teal_cmap():
    cdict = {
        'red': [(0.0, 1.0, 1.0), (0.2, 0.8, 0.8), (0.8, 0.3, 0.3), (1.0, 0.3, 0.3)],  # Beyazdan pastel teal'e geçiş
        'green': [(0.0, 1.0, 1.0), (0.2, 0.9, 0.9), (0.8, 0.7, 0.7), (1.0, 0.7, 0.7)],  # Pastel yeşil tonları
        'blue': [(0.0, 1.0, 1.0), (0.2, 0.9, 0.9), (0.8, 0.8, 0.8), (1.0, 0.8, 0.8)],  # Pastel mavi tonları
    }
    return mcolors.LinearSegmentedColormap('pastel_teal_cmap', segmentdata=cdict, N=256)


# Pastel turuncu renk paletini oluşturma
def create_pastel_orange_cmap(color):
    cdict = {
        'red': [(0.0, 1.0, 1.0), (0.2, 1.0, 1.0), (0.8, 1.0, 1.0), (1.0, color[0], color[0])],
        'green': [(0.0, 1.0, 1.0), (0.2, 0.8, 0.8), (0.8, 0.5, 0.5), (1.0, 0.4, 0.4)],
        'blue': [(0.0, 1.0, 1.0), (0.2, 0.5, 0.5), (0.8, 0.3, 0.3), (1.0, 0.1, 0.1)],
    }
    return mcolors.LinearSegmentedColormap('pastel_cmap', segmentdata=cdict, N=256)



# Veriyi uzun formata dönüştürme
def prepare_data_for_facetgrid(df):
    df_long = df.reset_index().melt(id_vars='Species', var_name='Condition', value_name='Percentage')
    return df_long



# FacetGrid kullanarak alt alta çizme
def plot_heatmaps(data_high, data_low, output, vmax_low):

    pastel_teal_cmap = create_pastel_teal_cmap()  # Teal tones
    pastel_orange_cmap = create_pastel_orange_cmap((1.0, 0.7, 0.3))  # Orange tones

    fig, axes = plt.subplots(nrows=2, figsize=(10, 12))  # 2 satırlı bir düzen
    
    # Yüksek değerler için heatmap
    sns.heatmap(data_high.pivot(index='Species', columns='Condition', values='Percentage'),
                cmap=pastel_orange_cmap, vmin=0, vmax=100, annot=True, fmt='.1f', linewidths=0.5, linecolor='gray',
                ax=axes[0], cbar_kws={'label': 'Percentage'})
    axes[0].set_title("Species Distribution Across Conditions (Percentage > 15)")
    
    # Düşük değerler için heatmap
    sns.heatmap(data_low.pivot(index='Species', columns='Condition', values='Percentage'),
                cmap=pastel_teal_cmap, vmin=0, vmax=vmax_low, annot=True, fmt='.1f', linewidths=0.5, linecolor='gray',
                ax=axes[1], cbar_kws={'label': 'Percentage'})
    axes[1].set_title("Species Distribution Across Conditions (Percentage <= 15)")
    
    plt.tight_layout()
    plt.savefig(output)  # Grafik dosyaya kaydedilir
    plt.show()



# Veriyi tek bir DataFrame'de birleştirme
def prepare_combined_data(df):
    # DataFrame'i uzun formata dönüştürme
    df_long = df.reset_index().melt(id_vars='Species', var_name='Condition', value_name='Percentage')
    return df_long


# Matrix formatında kaydetme
def save_combined_matrix(df, filename):
    # Veriyi matrix formatında kaydetme
    df.pivot_table(index='Species', columns='Condition', values='Percentage').to_csv(filename, sep='\t')



def main_heatmap(file_path, output ,combined_matrix_file_path):
    if not os.path.exists(file_path):
        print(f"Warning: File '{file_path}' does not exist.")
        return 0
    
    with open(file_path, 'r') as file:
        data = file.read()

    # Parse işlemi
    conditions = []
    species_data = []
    current_condition = None


    #Parsing .rcp file
    for line in data.splitlines():
        line = line.strip()
        if line.startswith(">"):
            current_condition = line[1:]
        elif line.startswith("?"):
            continue
        elif line == "---":
            current_condition = None
        elif current_condition:
            species = '_'.join((line.split(', ')[0],line.split(', ')[1]))
            percentage = float(line.split(', ')[2])
            print(species, percentage)
            conditions.append(current_condition)
            species_data.append((species, percentage))

    # DataFrame'e çevirme
    df = pd.DataFrame(species_data, columns=["Species", "Percentage"])
    df['Condition'] = conditions

    # Heatmap verisi oluşturma (pivot)
    heatmap_data = df.pivot(index="Species", columns="Condition", values="Percentage")

    # Eksik değerleri 0 ile doldurma
    heatmap_data = heatmap_data.fillna(0)

    # Verileri filtreleme
    heatmap_data_high = heatmap_data[heatmap_data.max(axis=1) > 15]
    heatmap_data_low = heatmap_data[heatmap_data.max(axis=1) <= 15]

    # Düşük yüzdelikler için vmax değerini belirleme
    vmax_low = heatmap_data_low.max().max()  # Düşük değerler için maksimum değer

   
    heatmap_data_high_long = prepare_data_for_facetgrid(heatmap_data_high)
    heatmap_data_low_long = prepare_data_for_facetgrid(heatmap_data_low)
    plot_heatmaps(heatmap_data_high_long, heatmap_data_low_long, output, vmax_low)
    combined_data_long = prepare_combined_data(heatmap_data)


    save_combined_matrix(combined_data_long, combined_matrix_file_path)


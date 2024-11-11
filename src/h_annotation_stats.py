import os
import pandas as pd
from Bio import SeqIO


def main_is_stat(data_dir, out_file):

    #current_dir = os.getcwd()
    current_dir = data_dir
    folders = [f for f in os.listdir(current_dir) if os.path.isdir(os.path.join(current_dir, f))]
    out_text_list = []
    out_text_segment = ""

    for folder in folders:

        the_condition = f"{folder}"
        #fasta_file = os.path.join(current_dir, folder, f"{folder}.fasta")
        tab_file = os.path.join(current_dir, folder, f"best_alignment.tab")

        if not os.path.exists(tab_file):
            print(f"Warning: '{tab_file}' not found. Skipping.")
            continue

        print('Processing the files')
        print(os.path.basename(tab_file))
        #print(os.path.basename(fasta_file))

        df_raw = pd.read_csv(tab_file, sep='\t')
        df = df_raw[df_raw['Note'] == 'IS_DB'].copy()
        df[['family_name', 'special_name']] = df['Subject ID'].str.split('_', expand=True)

        special_name_dict = {}

        for _, row in df.iterrows():
            special_name = row['special_name']
            family_name = row['family_name']
            special_name_dict[special_name] = family_name

        special_name_dict['no_trans'] = 'No_IS'
        print(special_name_dict)

        family_name_counts = df['family_name'].value_counts()
        special_name_counts = df['special_name'].value_counts()

        insert_count = df.shape[0]

        total_family_name_counts = family_name_counts.sum()
        total_special_name_counts = special_name_counts.sum()
        family_name_percentages = (family_name_counts / total_family_name_counts) * 100
        special_name_percentages = (special_name_counts / total_special_name_counts) * 100

        family_name_percentage_data = family_name_percentages.to_dict()
        special_name_percentage_data = special_name_percentages.to_dict()

        print("Family Name Frekansları ve Yüzdeleri:")
        for name, percentage in family_name_percentage_data.items():
            count = family_name_counts[name]  # Get the count for each family name
            print(f"{name}: {percentage:.2f}% ({count})")

        print("\nSpecial Name Frekansları ve Yüzdeleri:")

        out_text_segment = f">{the_condition}\n?{insert_count}\n"
        for name, percentage in special_name_percentage_data.items():
            count = special_name_counts[name]  # Get the count for each special name
            print(f"{name}: {percentage:.2f}% ({count})")
            out_text_segment = out_text_segment + f"{special_name_dict[name]}, {name}, {percentage:.2f}%, {count}\n"

        # Output file
        out_text_list.append(out_text_segment)

        print(f"\nInsert Count: {insert_count}")

    out_text = "---\n".join(out_text_list)

    with open(out_file, "w") as f:
        f.write(out_text)

#out_file = 'stats.rcp'
#main_is_stat(out_file)


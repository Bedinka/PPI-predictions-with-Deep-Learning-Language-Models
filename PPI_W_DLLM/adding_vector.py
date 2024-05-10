import pandas as pd
import os

df = pd.read_csv('your_file.tsv', sep='\t')

pickle_ca_dir = './Matrices_CA/'
pickle_mean_dir = './Matrices_Mean/'

for index, row in df.iterrows():
    protein_ids = row['Protein ID'].split('_')
    for protein_id in protein_ids:
        pickle_filename = f'm_ca_{protein_id}.pickle'
        pickle_path = os.path.join(pickle_ca_dir, pickle_filename)
        if os.path.exists(pickle_path):
            print(f'Pickle file found for protein ID {protein_id}: {pickle_path}')
        else:
            print(f'Pickle file not found for protein ID {protein_id}')
    vector_ca = f"{chains_CA[chainID1].dssp_index}_{chains_CA[chainID2].dssp_index}" 
# From HURI : http://www.interactome-atlas.org/download  HuRI.tsv
# Mapped in : https://www.uniprot.org/id-mapping from Ensembl To Uni/Swiss
# dowloaded as FASTA and interactions as TSV


from Bio import SeqIO

# Parse the FASTA file and extract protein identifiers
#
fasta_file = '/home/dina/Downloads/idmapping_2024_05_28 (1).fasta'
uniprot_ids = []

for record in SeqIO.parse(fasta_file, "fasta"):
    # Extract the UniProt identifier (e.g., 'J3KPH0', 'Q6NUM6', etc.)
    uniprot_id = record.id.split('|')[1]
    uniprot_ids.append(uniprot_id)

print(uniprot_ids)

import pandas as pd

# Load the HuRI dataset
huri_url = 'http://www.interactome-atlas.org/data/HuRI.tsv'
huri_df = pd.read_csv(huri_url, sep='\t', header=None, names=['Protein_A', 'Protein_B'])

# Load the mapping file from UniProt
mapping_file = '/home/dina/Downloads/idmapping_2024_05_28.tsv'
mapping_df = pd.read_csv(mapping_file, sep='\t')

uniport_column = 'Entry'  # Adjust this based on the actual column name for UniProt IDs in your mapping DataFrame

# Create a dictionary for quick lookup from identifiers in the 'From' column to UniProt ID
mapping_dict = pd.Series(mapping_df[uniport_column].values, index=mapping_df['From']).to_dict()


# Map Ensembl Gene IDs in HuRI to UniProt IDs
huri_df['Protein_A_UniProt'] = huri_df['Protein_A'].map(mapping_dict)
huri_df['Protein_B_UniProt'] = huri_df['Protein_B'].map(mapping_dict)

# Drop rows with missing mappings
huri_df.dropna(subset=['Protein_A_UniProt', 'Protein_B_UniProt'], inplace=True)

print(huri_df.head())

# Convert the list of UniProt IDs from the FASTA file to a set for faster lookups
alphafold_proteins = set(uniprot_ids)

# Filter interactions where both proteins are in the AlphaFold dataset
huri_df_filtered = huri_df[huri_df['Protein_A_UniProt'].isin(alphafold_proteins) & 
                           huri_df['Protein_B_UniProt'].isin(alphafold_proteins)]

# Save or print the filtered interactions
huri_df_filtered.to_csv('filtered_interactions.csv', index=False)
print(huri_df_filtered)
positive_pairs = huri_df_filtered[['Protein_A_UniProt', 'Protein_B_UniProt']].values.tolist()

with open('positive_pairs.txt', 'w') as f:
    for pair in positive_pairs:
        f.write(f"{pair[0]}\t{pair[1]}\n")
positive_pairs = huri_df_filtered[['Protein_A_UniProt', 'Protein_B_UniProt']].values.tolist()

import random

all_proteins = list(alphafold_proteins)
positive_pairs_set = set([tuple(sorted(pair)) for pair in positive_pairs])

negative_pairs = []

while len(negative_pairs) < len(positive_pairs):
    pair = random.sample(all_proteins, 2)
    if tuple(sorted(pair)) not in positive_pairs_set:
        negative_pairs.append(pair)

with open('negative_pairs.txt', 'w') as f:
    for pair in negative_pairs:
        f.write(f"{pair[0]}\t{pair[1]}\n")


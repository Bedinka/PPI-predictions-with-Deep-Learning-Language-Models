import torch
import numpy as np
import feature_extraction

# Call the main function of feature_extraction.py
feature_extraction.main()

# Access the protein_data dictionary

one_hot_sequences = []
distance_matrices_tensor = []
submatrices_tensor = []
rsa_values_tensor = []

# Convert protein sequences to one-hot encoding
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
num_amino_acids = len(amino_acids)

for protein_id, chain_data in feature_extraction.protein_data.items():
    for data in chain_data.items():
        sequences = data['residue_data']
        max_seq_length = max(len(seq) for seq in sequences)
        one_hot_seq = np.zeros((len(sequences), max_seq_length, num_amino_acids), dtype=np.float32)
        for i, seq in enumerate(sequences):
            for j, aa in enumerate(seq):
                if aa in amino_acids:
                    one_hot_seq[i, j, amino_acids.index(aa)] = 1
        one_hot_sequences.append(one_hot_seq)

        # Convert distance matrices to tensor
        distance_matrix = data['distance_matrix_CA']
        distance_matrices_tensor.append(torch.tensor(distance_matrix, dtype=torch.float32))

        # Convert submatrices to tensor
        submatrices_data = data['sub_matrixes']
        submatrices_tensor.append(torch.tensor(submatrices_data, dtype=torch.float32))

        # Convert RSA values to tensor
        rsa_values_data = data['residue_data']['rsa_value']
        rsa_values_tensor.append(torch.tensor(rsa_values_data, dtype=torch.float32))

# Concatenate lists into tensors
one_hot_sequences_tensor = np.concatenate(one_hot_sequences, axis=0)
distance_matrices_tensor = torch.cat(distance_matrices_tensor, dim=0)
submatrices_tensor = torch.cat(submatrices_tensor, dim=0)
rsa_values_tensor = torch.cat(rsa_values_tensor, dim=0)

# Print shapes for verification
print("One-hot encoded sequences shape:", one_hot_sequences_tensor.shape)
print("Distance matrices tensor shape:", distance_matrices_tensor.shape)
print("Submatrices tensor shape:", submatrices_tensor.shape)
print("RSA values tensor shape:", rsa_values_tensor.shape)
import feature_extraction
import numpy as np
import torch

# Call the main function of feature_extraction.py
feature_extraction.main()

# Extract the required data from feature_extraction.py or load from file if available
protein_sequences = feature_extraction.protein_data


# Convert protein sequences to one-hot encoding
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
num_amino_acids = len(amino_acids)
max_seq_length = max(len(seq) for seq in protein_sequences)
one_hot_sequences = np.zeros((len(protein_sequences), max_seq_length, num_amino_acids), dtype=np.float32)
for i, seq in enumerate(protein_sequences):
    for j, aa in enumerate(seq):
        if aa in amino_acids:
            one_hot_sequences[i, j, amino_acids.index(aa)] = 1

# Convert distance matrices, submatrices, and RSA values to tensors
distance_matrices_tensor = torch.tensor(distance_matrices, dtype=torch.float32)
submatrices_tensor = torch.tensor(submatrices, dtype=torch.float32)
rsa_values_tensor = torch.tensor(rsa_values, dtype=torch.float32)

# Print shapes for verification
print("One-hot encoded sequences shape:", one_hot_sequences.shape)
print("Distance matrices tensor shape:", distance_matrices_tensor.shape)
print("Submatrices tensor shape:", submatrices_tensor.shape)
print("RSA values tensor shape:", rsa_values_tensor.shape)

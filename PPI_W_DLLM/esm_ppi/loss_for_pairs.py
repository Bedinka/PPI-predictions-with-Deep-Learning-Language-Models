import numpy as np
from scipy.optimize import linear_sum_assignment
from transformers import AutoTokenizer, EsmForMaskedLMpython 

import torch
import csv

tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
model = EsmForMaskedLM.from_pretrained("facebook/esm2_tesm2_t6_8M_UR50D6_8M_UR50D")

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = model.to(device)

all_proteins = []
tsv_file_path = "protein_sequences.tsv"
for i in range (6):
    with open(tsv_file_path, "r", newline="") as tsv_file:
        # Create a CSV reader object
        tsv_reader = csv.reader(tsv_file, delimiter="\t")
        
        # Iterate over each row in the TSV file
        for row in tsv_reader:
            # Assuming the protein sequence is in the first column (index 0) of each row
            protein_sequence = row[0]
            all_proteins.append(protein_sequence)

BATCH_SIZE = 2
NUM_MASKS = 10
P_MASK = 0.15

# Function to compute MLM loss for a batch of protein pairs
def compute_mlm_loss_batch(pairs):
    avg_losses = []
    for _ in range(NUM_MASKS):
        # Tokenize the concatenated protein pairs
        inputs = tokenizer(pairs, return_tensors="pt", truncation=True, padding=True, max_length=1022)
        
        # Move input tensors to GPU if available
        inputs = {k: v.to(device) for k, v in inputs.items()}
        
        # Get the mask token ID
        mask_token_id = tokenizer.mask_token_id
        
        # Clone input IDs for labels
        labels = inputs["input_ids"].clone()
        
        # Randomly mask 15% of the residues for each sequence in the batch
        for idx in range(inputs["input_ids"].shape[0]):
            mask_indices = np.random.choice(inputs["input_ids"].shape[1], size=int(P_MASK * inputs["input_ids"].shape[1]), replace=False)
            inputs["input_ids"][idx, mask_indices] = mask_token_id
            labels[idx, [i for i in range(inputs["input_ids"].shape[1]) if i not in mask_indices]] = -100
        
        # Compute the MLM loss
        outputs = model(**inputs, labels=labels)
        avg_losses.append(outputs.loss.item())
    
    # Return the average loss for the batch
    return sum(avg_losses) / NUM_MASKS

# Compute loss matrix
loss_matrix = np.zeros((len(all_proteins), len(all_proteins)))

for i in range(len(all_proteins)):
    for j in range(i+1, len(all_proteins), BATCH_SIZE):  # to avoid self-pairing and use batches
        pairs = [all_proteins[i] + all_proteins[k] for k in range(j, min(j+BATCH_SIZE, len(all_proteins)))]
        batch_loss = compute_mlm_loss_batch(pairs)
        for k in range(len(pairs)):
            loss_matrix[i, j+k] = batch_loss
            loss_matrix[j+k, i] = batch_loss  # the matrix is symmetric

# Set the diagonal of the loss matrix to a large value to prevent self-pairings
np.fill_diagonal(loss_matrix, np.inf)

rows, cols = linear_sum_assignment(loss_matrix)
optimal_pairs = list(zip(rows, cols))

print(optimal_pairs)
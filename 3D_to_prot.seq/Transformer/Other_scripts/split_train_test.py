import random

def remove_random_sequences(file1_path, file2_path, removed_file1_path, removed_file2_path, percentage):
    sequences1 = []
    sequences2 = []

    with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
        sequences1 = file1.readlines()
        sequences2 = file2.readlines()

    num_sequences = min(len(sequences1), len(sequences2))
    num_removed_sequences = int(num_sequences * percentage)

    random.seed(42)  # Set a fixed seed for reproducibility
    removed_indices = random.sample(range(num_sequences), num_removed_sequences)
    removed_indices.sort(reverse=True) 

    with open(file1_path, 'w') as file1, open(file2_path, 'w') as file2, \
            open(removed_file1_path, 'w') as removed_file1, open(removed_file2_path, 'w') as removed_file2:
        for i in range(num_sequences):
            if i in removed_indices:
                removed_file1.write(sequences1[i])
                removed_file2.write(sequences2[i])
            else:
                file1.write(sequences1[i])
                file2.write(sequences2[i])

file1_path = './data/train_fasta.90.txt'
file2_path = './data/train_structure.90.txt'
removed_file1_path = './data/Test/test_fasta.txt'
removed_file2_path = './data/Test/test_structure.txt'
percentage = 0.1

remove_random_sequences(file1_path, file2_path, removed_file1_path, removed_file2_path, percentage)

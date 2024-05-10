# Read the representative protein sequences from the second file
representative_sequences = []
with open('interactom_nonredundant', 'r') as file:
    lines = file.readlines()
    sequence = ''
    for line in lines:
        if line.startswith("['") and '_' in line:
            sequence = line.strip()[2:-3]  # Extract sequence from the line
            representative_sequences.append(sequence)

print("Representative sequences:", representative_sequences)

# Read the first file and filter lines based on residue sequences in the fourth column
filtered_lines = []
with open('bert_train_5 (copy).tsv', 'r') as file:
    lines = file.readlines()
    for line in lines:
        parts = line.strip().split('\t')
        if len(parts) >= 4:
            residue_sequence = parts[3]
            if residue_sequence in representative_sequences:
                filtered_lines.append(line)
        else:
            print(f"Issue with line: {line}")  # Print lines that do not meet the expected format

print("Filtered lines:", filtered_lines)

# Write the filtered lines to a new file
with open('filtered_first_file.tsv', 'w') as file:
    file.writelines(filtered_lines)

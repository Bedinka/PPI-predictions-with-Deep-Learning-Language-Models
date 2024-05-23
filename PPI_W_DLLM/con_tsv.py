import pandas as pd

# Paths to the TSV files
file1 = 'concatenated_file.tsv'
file2 = 'bert_train_with_vector_002003.tsv'

df1 = pd.read_csv(file1, sep='\t')
df2 = pd.read_csv(file2, sep='\t')

# Concatenate the dataframes, ignoring the header row in the second dataframe
concatenated_df = pd.concat([df1, df2], ignore_index=True)

# Save the concatenated dataframe to a new TSV file
output_file = 'concatenated_file_v1.tsv'
concatenated_df.to_csv(output_file, sep='\t', index=False)

print(f"Files {file1} and {file2} concatenated and saved to {output_file}")
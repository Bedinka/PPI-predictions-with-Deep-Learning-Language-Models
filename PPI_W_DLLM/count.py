import pandas as pd

def count_ones_and_zeros(file_path):
    try:
        # Read the TSV file using pandas
        df = pd.read_csv(file_path, sep='\t')
        
        # Display the first few rows to ensure the file is read correctly
        print("First few rows of the dataframe:")
        print(df.head())

        # Ensure the column exists
        if 'Interact' not in df.columns:
            print("The column 'Interact' does not exist in the file.")
            return

        # Count the occurrences of '1' and '0' in the 'Interact' column
        count_ones = (df['Interact'] == 1).sum()
        count_zeros = (df['Interact'] == 0).sum()

        print(f'Number of 1s: {count_ones}')
        print(f'Number of 0s: {count_zeros}')
    except Exception as e:
        print(f"An error occurred: {e}")

# Path to your CSV file
file_path = '2024-06-04_17-36-41/bert_train_with_vector_2024-06-04_17-36-41_500_s_wDataloader.tsv'

# Get the counts
count_ones_and_zeros(file_path)

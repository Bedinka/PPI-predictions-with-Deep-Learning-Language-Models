import pandas as pd

def read_tsv_file(file_path):
    try:
        # Read the TSV file using pandas
        df = pd.read_csv(file_path, sep='\t', engine='python')
        return df
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None

def show_column_info(dataset_name, dataset):
    print(f"Data types for the {dataset_name} dataset:")
    print(dataset.info())
    print('===' * 14)

def show_dataset_summary_statistics(datasets):
    for name, data in datasets.items():
        print(f"{name.capitalize()} dataset summary statistics:")
        print('---' * 15)
        print(data.describe())
        print('===' * 20)
        print()

def show_missing_values(datasets):
    for name, data in datasets.items():
        print(f"Missing values in the {name.capitalize()} dataset:")
        print(data.isnull().sum())
        print('===' * 18)
        print()

# Specify the path to your TSV file here
file_path = '/home/dina/Documents/PPI_WDLLM/2024-06-05_15-23-05/bert_train_with_vector_2024-06-05_15-23-05_300_s_wDataloader_negativex10_fortesting.tsv'

# Read the dataset
df = read_tsv_file(file_path)

if df is not None:
    print(f"Test Dataset: {df.shape}")
    show_column_info('Test', df)
    print()

    # Checking for the summary statistics of the datasets
    datasets = {'test': df}
    show_dataset_summary_statistics(datasets)

    # Check for missing values in the datasets
    show_missing_values(datasets)
else:
    print("Failed to read the dataset.")


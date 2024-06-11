import pandas as pd
import os
import numpy as np
import pickle
import test_autoencoder

def save_vectors_as_pickle(pickle_dir, sub_size, input_tsv_for_vec, auto_model_dir, encode_test_model_name, vector_pickle_dir):
    df = pd.read_csv(input_tsv_for_vec, sep='\t')
    auto_model_path = os.path.join(auto_model_dir, encode_test_model_name)
    
    for index, row in df.iterrows():
        protein_ids = row['Protein ID'].split('_')
        
        pickle_1 = f'm_ca_{protein_ids[0]}.pickle'
        path = os.path.join(pickle_dir, pickle_1)
        encoded_vector_1 = test_autoencoder.main(auto_model_path, sub_size, path)
        
        pickle_2 = f'm_ca_{protein_ids[1]}.pickle'
        path = os.path.join(pickle_dir, pickle_2)
        encoded_vector_2 = test_autoencoder.main(auto_model_path, sub_size, path)
        
        if encoded_vector_1 is not None and encoded_vector_2 is not None:
            vector = np.concatenate((encoded_vector_1, encoded_vector_2))
            vector_pickle_path = os.path.join(vector_pickle_dir, f'vector_{index}.pickle')
            with open(vector_pickle_path, 'wb') as f:
                pickle.dump(vector, f)

if __name__ == "__main__":
    pickle_dir = '/home/dina/Documents/PPI_WDLLM/Matrices_CA/train'
    sub_size = 7 
    input_tsv_for_vec = '2024-06-05_15-23-05/filtered_file_train_2024-06-05_15-23-05_300_s_wDataloader_negativex10_fortesting.tsv'
    encode_test_model_name = '/home/dina/Documents/PPI_WDLLM/autoencoder_dina_models/autoencoder_trained_2024-06-04_17-36-41_500_s_wDataloader_2.keras'
    vector_pickle_dir = '/home/dina/Documents/PPI_WDLLM/Vector_Pickles'
    os.makedirs(vector_pickle_dir, exist_ok=True)
    save_vectors_as_pickle(pickle_dir, sub_size, input_tsv_for_vec, auto_model_dir, encode_test_model_name, vector_pickle_dir)

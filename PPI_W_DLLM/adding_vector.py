import pandas as pd
import os
import numpy as np
import test_autoencoder


def main(pickle_dir, sub_size, input_tsv_for_vec, auto_model_dir,encode_test_model_name, vectorized_tsv ):

    vector_ca_list = []
    df = pd.read_csv(input_tsv_for_vec, sep='\t')
    auto_model_path = os.path.join(auto_model_dir, encode_test_model_name)
    for index, row in df.iterrows():
        protein_ids = row['Protein ID'].split('_')
        pickle_1 = f'm_ca_{protein_ids[0]}.pickle'
        path = os.path.join(pickle_dir, pickle_1)
        encoded_vector_1= test_autoencoder.main(auto_model_path, sub_size, path)

        pickle_2 = f'm_ca_{protein_ids[1]}.pickle'
        path = os.path.join(pickle_dir, pickle_2)
        encoded_vector_2= test_autoencoder.main(auto_model_path, sub_size, path)

        ca1, ca2 = None, None
        if encoded_vector_1 is not None:
            print(f'Pickle file 1 found for protein ID {protein_ids[0]}')
            ca1 = np.array(encoded_vector_1)
        else:
            print(f'Pickle file not found for protein ID {protein_ids[0]}')
        if encoded_vector_2 is not None:
            print(f'Pickle file 2 found for protein ID {protein_ids[1]}')
            ca2 = np.array([encoded_vector_1])
        else:
            print(f'Pickle file not found for protein ID {protein_ids[1]}')

        if ca1 is not None and ca2 is not None:
            vector_ca = f"{ca1}_{ca2}"
        else:
            vector_ca = None  
        vector_ca_list.append(vector_ca)

    df['Vector'] = vector_ca_list
    df.to_csv(vectorized_tsv, sep='\t', index=False)

if __name__ == "__main__":
    pickle_dir = ''
    sub_size = 7 
    path = ''
    encode_test_model_name = ''
    vectorized_tsv = 'bert_test_with_vector_005.tsv'
    main(pickle_dir, sub_size, path, auto_model_dir, encode_test_model_name, vectorized_tsv )

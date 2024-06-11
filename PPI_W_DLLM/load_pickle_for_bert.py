import pandas as pd
import os
import pickle

def load_vectors_for_bert(input_tsv, vector_pickle_dir):
    df = pd.read_csv(input_tsv, sep='\t')
    vectors = []

    for index, row in df.iterrows():
        vector_pickle_path = os.path.join(vector_pickle_dir, f'vector_{index}.pickle')
        if os.path.exists(vector_pickle_path):
            with open(vector_pickle_path, 'rb') as f:
                vector = pickle.load(f)
            vectors.append(vector)
        else:
            vectors.append(None)  # Handle missing vectors appropriately
    
    df['Vector'] = vectors
    return df

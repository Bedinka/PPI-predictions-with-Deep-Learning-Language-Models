import pandas as pd
import os
import numpy as np
import data_autoencoding

df = pd.read_csv('bert_big_test.tsv', sep='\t')
pickle_ca_dir = './Matrices_CA_test/'
pickle_mean_dir = './Matrices_Mean_test/'
size = 7 
vector_ca_list = []

model_name = 'dina_model_final.keras'
processed_sample = 2
latent_dim = 2
size = 7
epoch = 10
i = 0
TEST =True
SAVE=True
batch_size =  32

for index, row in df.iterrows():
    protein_ids = row['Protein ID'].split('_')
    pickle_1 = f'm_ca_{protein_ids[0]}.pickle'
    unseen_data_path = os.path.join(pickle_ca_dir, pickle_1)
    encoded_vector_1= data_autoencoding.main(latent_dim, model_name, processed_sample , size , SAVE , epoch , TEST, unseen_data_path, batch_size )

    pickle_2 = f'm_ca_{protein_ids[1]}.pickle'
    unseen_data_path = os.path.join(pickle_ca_dir, pickle_2)
    encoded_vector_2= data_autoencoding.main(latent_dim, model_name, processed_sample , size , SAVE , epoch , TEST, unseen_data_path, batch_size )

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
df.to_csv('bert_train_with_vector_002003.tsv', sep='\t', index=False)

import os
import sys
import datetime
import logging
import pandas as pd
import warnings
from Bio.PDB.PDBParser import PDBConstructionWarning

import setup_run

logging.basicConfig(filename='script_log.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
# setting what script to run 
F1_RUN = False
F2_RUN = True 
S_REMOVE = True
PDB2FASTA = True
A_RUN = True
A_TEST = True
VECTOR = False
VECTOR_SAVE = True
B_RUN = True
B_TEST = False

# Set up step information
steps = [
    ("Feature extraction for single chain", F1_RUN),
    ("Feature extraction for double chain", F2_RUN),
    ("PDB to fasta conversion", PDB2FASTA),
    ("Removing redundant sequences", S_REMOVE),
    ("Autoencoder training", A_RUN),
    ("Autoencoder testing", A_TEST),
    ("Vector append", VECTOR),
    ("Vector save", VECTOR_SAVE),
    ("BERT training", B_RUN),
    ("BERT testing", B_TEST)
]

#Extra valuese for use 
time = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
info = '1300_s_wDataloader_interactom'

set = True # if true it will set up a new run directory 
path = ''
# Setup the run directory
    # Setup LoggerWriter to redirect print statements to logging
run_directory = setup_run.setup_run_directory(set, path )

#PDB files directory
pdb_file_dir = 'TRAIN'
logging.info("Trainin directory: %s ", pdb_file_dir)
#Autoencoder models directory
auto_model_dir = '/home/dina/Documents/PPI_WDLLM/autoencoder_dina_models'

# Setting pickle directory
pickle_dir_ca = '/home/dina/Documents/PPI_WDLLM/Matrices_CA/'
pickle_dir_mean = '/home/dina/Documents/PPI_WDLLM/Matrices_Mean/'
p_file_dir_train = 'train'
p_file_dir_test = 'test'

# CHOOSE CA OR MEAN for amino acid  for computing 
pickle_dir= pickle_dir_ca
ca_dir_path = os.path.join(pickle_dir_ca, p_file_dir_train)
mean_dir_path = os.path.join(pickle_dir_mean, p_file_dir_train)
pickle_train_path = os.path.join(pickle_dir, p_file_dir_train)
pickle_test_path = os.path.join(pickle_dir, p_file_dir_test)

logging.info("######################" * 10) 
logging.info("Run TRAIN with pickle directory : %s" , pickle_train_path)
logging.info("######################"* 10) 
logging.info("Run TEST with pickle directory : %s" , pickle_test_path)

for step_name, should_run in steps:
    if should_run:
        logging.info("Starting step: %s", step_name)
        print("Starting step:", step_name)
        try:
            if step_name == "Feature extraction for single chain":
                positive_pairs_txt = '/home/dina/Documents/PPI_WDLLM/positive_pairs.txt' #% (time, info)
                if F1_RUN:
                    import feature_extraction_v2
                    logging.info("-" * 100)
                    logging.info("Starting feature extraction for SINGLE chain with:\n"
                        "sample_batch=%d\n"
                        "sub_size=%d\n"
                        "feature_tsv=%s", sample, sub_size, feature_tsv_output_name)
                    feature_extraction_v2.main(sample, sub_size, feature_tsv_output_name, ca_dir_path, mean_dir_path, pdb_file_dir, positive_pairs_txt)
                pass
            elif step_name == "Feature extraction for double chain":
                if F2_RUN:
                import feature_extraction
                logging.info("-" * 100)
                logging.info("Starting feature extraction for DOUBLE chain with:\n"
                    "sample_batch=%d\n"
                    "sub_size=%d\n"
                    "feature_tsv=%s", sample, sub_size, feature_tsv_output_name)
                feature_extraction.main(sample, sub_size, feature_tsv_output_name, ca_dir_path, mean_dir_path, pdb_file_dir)

            elif step_name == "PDB to fasta conversion":
                # Insert code for PDB to fasta conversion
                pass
            elif step_name == "Removing redundant sequences":
                # Insert code for removing redundant sequences
                pass
            elif step_name == "Autoencoder training":
                # Insert code for autoencoder training
                pass
            elif step_name == "Autoencoder testing":
                # Insert code for autoencoder testing
                pass
            elif step_name == "Vector append":
                # Insert code for vector append
                pass
            elif step_name == "Vector save":
                # Insert code for vector save
                pass
            elif step_name == "BERT training":
                # Insert code for BERT training
                pass
            elif step_name == "BERT testing":
                # Insert code for BERT testing
                pass
            logging.info("Step %s completed successfully", step_name)
            print("Step completed successfully")
        except Exception as e:
            logging.error("Error occurred in step %s: %s", step_name, str(e))
            print("Error occurred in step:", step_name)
            break 
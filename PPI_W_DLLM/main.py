import os
import sys
import datetime
import logging
import pandas as pd

import feature_extraction
import train_autoencoding
import test_autoencoder
import combine_input_bert
import redundancy_remove
import adding_vector
import pdb2fatsa
import cd_hit

import warnings
from Bio.PDB.PDBParser import PDBConstructionWarning

# Set up a filter to ignore warnings from the pdb2fasta module
warnings.filterwarnings("ignore", category=PDBConstructionWarning)

# Setup LoggerWriter to redirect print statements to logging
class LoggerWriter:
    def __init__(self, level):
        self.level = level

    def write(self, message):
        if message.strip():
            self.level(message)

    def flush(self):
        pass

# Initialize logging
run_directory = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
if not os.path.exists(run_directory):
    os.makedirs(run_directory)
    

log_file = os.path.join(run_directory, 'run.log')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', handlers=[
    logging.FileHandler(log_file),
    logging.StreamHandler()
])

sys.stdout = LoggerWriter(logging.info)
sys.stderr = LoggerWriter(logging.error)

logging.info("Run directory created at %s", run_directory)
#os.chdir(run_directory)

F_RUN = False 
S_REMOVE = False
PDB2FASTA = False
A_RUN = True
A_TEST = False
VECTOR = False
B_RUN = False
B_TEST = False

pdb = 'TRAIN'
auto_dir = './autoencoder_dina_models/'

# Setting pickle directory
pickle_dir_ca = '/home/dina/Documents/PPI_WDLLM/Matrices_CA/'
pickle_dir_mean = '/home/dina/Documents/PPI_WDLLM/Matrices_Mean/'
p_file = 'train'
pickle_ca_path = os.path.join(pickle_dir_ca, p_file)
pickle_mean_path = os.path.join(pickle_dir_mean, p_file)

# CHOOSE CA OR MEAN
pickle_dir = pickle_ca_path

# COMMON ATTRIBUTES

###############################################
# Running Feature extraction
sample = 2000
sub_size = 7
feature_tsv = "bert_train_2000_ca.tsv"

if F_RUN:
    logging.info("Starting feature extraction with sample=%d, sub_size=%d, feature_tsv=%s", sample, sub_size, feature_tsv)
    feature_extraction.main(sample, sub_size, feature_tsv, pickle_ca_path, pickle_mean_path, pdb)

###############################################
# PDB to fasta
if PDB2FASTA:
    interactom_fasta = 'interactom_train_2000_ca.fasta'
    logging.info("Running PDB to fasta conversion with pdb=%s, interactom_fasta=%s", pdb, interactom_fasta)
    pdb2fatsa.main(pdb, interactom_fasta)

    with open('/dev/null', 'w') as devnull:
        sys.stderr = devnull
        pdb2fatsa.main(pdb, interactom_fasta)
        sys.stderr = LoggerWriter(logging.error)

    output_nr_file = 'nonredundant'
    cd_hit.main(interactom_fasta , output_nr_file )


###############################################
# Removing non-representative sequences
if S_REMOVE:
    input_file = feature_tsv
    interactome_file = output_nr_file
    re_output_file = 'filtered_file_train_2000_ca.tsv'
    logging.info("Removing non-representative sequences with input_file=%s, interactome_file=%s, re_output_file=%s", input_file, interactome_file, re_output_file)
    redundancy_remove.main(input_file, interactome_file, re_output_file, run_directory) 

################################################
# Autoencoder TRAINing
# with multiple input attributes

# Train attributes
SAVE = True
batch_size = 32
pd_results = []
processed_sample_values = [2000]
size_values = [7]
latent_dim_values = [2]
epochs = [10]
ranges = 3

if A_RUN:
    logging.info("Starting Autoencoder training")
    for processed_sample in processed_sample_values:
        for latent_dim in latent_dim_values:
            for size in size_values:
                for epoch in epochs:
                    for i in range(ranges):
                        model_name = 'dina_model_sample_%d_dim_%d_size_%d_epochs_%d_index_%d.keras' % (processed_sample, latent_dim, size, epoch, i)
                        logging.info("Training model %s with latent_dim=%d, epoch=%d, processed_sample_num=%d, matrix_size=%dx%d", model_name, latent_dim, epoch, processed_sample, size, size)
                        try:
                            encoded_vector, collected_data = train_autoencoding.main(latent_dim, model_name, processed_sample, size, SAVE, epoch, batch_size, auto_dir, pickle_dir)
                            pd_results.append({
                                "Model": model_name,
                                "Latent Dimension": latent_dim,
                                "Matrix Size": size,
                                "Processed Sample": processed_sample,
                                "Epochs": epoch,
                                "Spearman Correlation": collected_data["spearman_correlation"],
                                "Spearman p-value": collected_data["spearman_p_value"],
                                "Spearman Correlation Dim1_2": collected_data["spearman_correlation_dim1_2"],
                                "Spearman p-value Dim1_2": collected_data["spearman_p_value_dim1_2"]
                            })
                        except IndexError as e:
                            logging.error("An IndexError occurred: %s. Skipping to the next iteration.", e)
                            continue
    df = pd.DataFrame(pd_results)
    tsv = 'autoencoder_train_sample_%d_dim_%d_size_%d_epochs_%d.tsv' % (processed_sample, latent_dim, size, epoch)
    df.to_csv(tsv, sep='\t', index=False)
    logging.info("Autoencoder training results saved to %s", tsv)

################################################
# Autoencoder TESTing

# Test attributes
if A_TEST:
    logging.info("Starting Autoencoder testing")
    pickle_data = [os.path.join(pickle_dir, file) for file in os.listdir(pickle_dir) if file.endswith('.pickle')]
    for p in pickle_data:
        encode_test_model_name = ''  # Define model name accordingly
        logging.info("Testing model %s with unseen data from %s", encode_test_model_name, p)
        try:
            encoded_vector_test = test_autoencoder.main(encode_test_model_name, sub_size, pickle_dir)
        except IndexError as e:
            logging.error("An IndexError occurred: %s. Skipping to the next iteration.", e)

################################################
# Vector append
if VECTOR:
    input_tsv_for_vec = ''  # Define input TSV file for vector
    vectorized_tsv_name = 'bert_test_with_vector_005'
    logging.info("Appending vector with input_tsv_for_vec=%s, encode_test_model_name=%s, vectorized_tsv_name=%s", input_tsv_for_vec, encode_test_model_name, vectorized_tsv_name)
    adding_vector.main(pickle_dir, sub_size, input_tsv_for_vec, encode_test_model_name, vectorized_tsv_name)

################################################
# BERT Training

output_stat_tsv = 'BERT_training_stats_withvectors_004.tsv'
train_input_tsv = 'concatenated_file_v1.tsv'  # Define vectorized TSV name

if B_RUN:
    info = 'rep01020304'
    combined_fields = ["Residues", "DSSP Structure", "DSSP Index", "Vector"]
    try:
        full_df = pd.read_csv(output_stat_tsv, sep='\t')
    except FileNotFoundError:
        full_df = pd.DataFrame()

    for i in range(1, len(combined_fields) + 1):
        bert_model_name = 'bert_attrnum_%d_no_%s' % (i, info)
        logging.info("Running BERT model %s with combined fields: %s", bert_model_name, combined_fields[:i])
        fields_to_combine = combined_fields[:i]
        df_stats = combine_input_bert.main(bert_model_name, train_input_tsv, fields_to_combine)
        df_stats.loc[0, 'model_name'] = bert_model_name
        full_df = pd.concat([full_df, df_stats], ignore_index=True)
    full_df.to_csv(output_stat_tsv, sep='\t', index=False)
    logging.info("BERT training stats saved to %s", output_stat_tsv)

################################################
# BERT Testing
if B_TEST:
    logging.info("Starting BERT testing")

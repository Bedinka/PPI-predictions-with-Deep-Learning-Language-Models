import os
import sys
import datetime
import logging
import pandas as pd
import warnings
from Bio.PDB.PDBParser import PDBConstructionWarning

import setup_run

# setting what script to run 
F1_RUN = False
F2_RUN = True 
S_REMOVE = True
PDB2FASTA = True
A_RUN = True
A_TEST = True
VECTOR = True
B_RUN = True
B_TEST = False

#Extra valuese for use 
time = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
info = 'pipelin_testing_alphafoldpdbs'
# Set up a filter to ignore warnings from the pdb2fasta module
warnings.filterwarnings("ignore", category=PDBConstructionWarning)
# Suppress TensorFlow logging
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
tf_logger = logging.getLogger('tensorflow')
tf_logger.setLevel(logging.ERROR)

# Setup the run directory
    # Setup LoggerWriter to redirect print statements to logging
run_directory = setup_run.setup_run_directory()

#PDB files directory
pdb_file_dir = 'TRAIN'
#Autoencoder models directory
auto_model_dir = 'autoencoder_dina_models'

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

logging.info("######################") 
logging.info("Run TRAIN with pickle directory : %s" , pickle_train_path)
logging.info("######################") 
logging.info("Run TEST with pickle directory : %s" , pickle_test_path)


##############################################################################################
#Feature extractions attributes: process sample number , sub matrix size (x,x)
sample = 100
sub_size = 7
feature_tsv_output_name  = "%s_%s.tsv" % (time, info)


###############################################
# Running Feature extraction for single pdb (1 chain), need a positiv pairs txt to create the interactions 


positive_pairs_txt = '/home/dina/Documents/PPI_WDLLM/positive_pairs.txt' #% (time, info)
if F1_RUN:
    import feature_extraction_v2
    logging.info("Starting feature extraction for SINGLE chain with:\n"
        "sample_batch=%d\n"
        "sub_size=%d\n"
        "feature_tsv=%s", sample, sub_size, feature_tsv_output_name)
    feature_extraction_v2.main(sample, sub_size, feature_tsv_output_name, ca_dir_path, mean_dir_path, pdb_file_dir, positive_pairs_txt)

###############################################
# Running Feature extraction for interactom pdb ( 2 chains)

if F2_RUN:
    import feature_extraction
    logging.info("Starting feature extraction for DOUBLE chain with:\n"
        "sample_batch=%d\n"
        "sub_size=%d\n"
        "feature_tsv=%s", sample, sub_size, feature_tsv_output_name)
    feature_extraction.main(sample, sub_size, feature_tsv_output_name, ca_dir_path, mean_dir_path, pdb_file_dir)

###############################################
# PDB to fasta
if PDB2FASTA:
    import pdb2fatsa
    import cd_hit
    input_tsv_p2f = feature_tsv_output_name
    interactom_fasta_output_name = 'interactom_train_%s_%s_ca.fasta' % (time, info)
    logging.info("Running PDB to fasta conversion with pdb=%s\n"", interactom_fasta=%s", pdb_file_dir, interactom_fasta_output_name)
    pdb2fatsa.main(input_tsv_p2f, pdb_file_dir, interactom_fasta_output_name)

    with open('/dev/null', 'w') as devnull:
        sys.stderr = devnull
        pdb2fatsa.main(input_tsv_p2f,pdb_file_dir, interactom_fasta_output_name)

    output_nr_file_name  = 'nonredundant_%s_%s' % (time, info)
    cd_hit.main(interactom_fasta_output_name , output_nr_file_name )


###############################################
# Removing non-representative sequences
if S_REMOVE:
    import redundancy_remove
    input_file = feature_tsv_output_name
    interactome_file = output_nr_file_name
    re_output_file = 'filtered_file_train__%s_%s.tsv'  % (time, info)
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
    import train_autoencoding
    logging.info("Starting Autoencoder training")
    for processed_sample in processed_sample_values:
        for latent_dim in latent_dim_values:
            for size in size_values:
                for epoch in epochs:
                    for i in range(ranges):
                        model_name = 'autoencoder_train__%s_%s_%d.keras' % (time, info, i) #'dina_model_sample_%d_dim_%d_size_%d_epochs_%d_index_%d.keras' % (processed_sample, latent_dim, size, epoch, i)
                        logging.info("Training model %s with latent_dim=%d,\n"" epoch=%d,\n"" processed_sample_num=%d,\n"" matrix_size=%dx%d", model_name, latent_dim, epoch, processed_sample, size, size)
                        try:
                            encoded_vector, collected_data = train_autoencoding.main(latent_dim, model_name, processed_sample, size, SAVE, epoch, batch_size, auto_model_dir, pickle_dir)
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
                            logging.info("Epoch %d completed for model %s", epoch, model_name)
                            logging.info("Autoencoder training results: %s", pd_results)
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
    import test_autoencoder
    logging.info("Starting Autoencoder testing")
    pickle_data = [os.path.join(pickle_dir, file) for file in os.listdir(pickle_dir) if file.endswith('.pickle')]
    best_model_name = ''
    correlation_best = 0
    for p in pickle_data[:3]:
        for i in range(ranges):
            encode_test_model_name = 'autoencoder_train_%s_%s_%d.keras' % (time, info, i)  # Define model name accordingly
            auto_model_path = os.path.join(auto_model_dir, encode_test_model_name)
            logging.info("Testing model %s with unseen data from %s", encode_test_model_name, p)
            try:
                encoded_vector_test, correlation_unseen = test_autoencoder.main(auto_model_path, sub_size, pickle_dir)
            except IndexError as e:
                logging.error("An IndexError occurred: %s. Skipping to the next iteration.", e)
        if correlation_unseen >= correlation_best:
            best_model_name = encode_test_model_name
            logging.info("Best model %s for %s file with correlation: %d ", encode_test_model_name, p,  correlation_best)

        

################################################
# Vector append
if VECTOR:
    import adding_vector
    input_tsv_for_vec = re_output_file  # Define input TSV file for vector
    vectorized_tsv_name = 'bert_train_with_vector_%s_%s.tsv' % (time, info)
    logging.info("Appending vector with input_tsv_for_vec=%s,\n"" encode_test_model_name=%s,\n"" vectorized_tsv_name=%s", input_tsv_for_vec, encode_test_model_name, vectorized_tsv_name)
    adding_vector.main(pickle_dir, sub_size, input_tsv_for_vec, auto_model_dir, best_model_name, vectorized_tsv_name)

################################################
# BERT Training

output_stat_tsv = 'BERT_training_stats_withvectors_%s_%s.tsv' % (time, info)
train_input_tsv = vectorized_tsv_name #'vectorized_tsv_name.tsv'  # Define vectorized TSV name

if B_RUN:
    import combine_input_bert
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

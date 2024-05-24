import feature_extraction
import data_autoencoding
import run_autoencoder
import combine_input_bert
import redundancy_remove
import adding_vector

import pandas as pd
import os
import datetime 

F_RUN = False 
B_RUN = False
B_TEST = False
S_REMOVE = False
A_RUN = True
A_TEST = False

# Setting run directory 
run_directory = datetime.datetime.now().strftime()
if not os.path.exists(run_directory):
        os.makedirs(run_directory)
        os.chdir(run_directory)
auto_dir = './autoencoder_dina_models/'

# Setting pickle directory
pickle_dir_ca = '/home/dina/Documents/PPI_WDLLM/Matrices_CA/'
pickle_dir_mean = '/home/dina/Documents/PPI_WDLLM/Matrices_Mean/'
p_file = 'test'
pickle_ca_path = os.path.join(pickle_dir_ca, p_file)
pickle_mean_path = os.path.join(pickle_dir_mean, p_file)

#CHOOSE CA OR MEAN
pickle_dir = pickle_ca_path

# COMMON ATTRIBUTES


###############################################
# Running Feature extraction
sample = 50
sub_size = 7
feature_tsv = "bert_small_test.tsv"

if F_RUN == True:
    feature_extraction.main(sample, sub_size, feature_tsv , pickle_ca_path, pickle_mean_path)


###############################################
# Removing non representative sequences

if S_REMOVE == True :
    input_file = feature_tsv
    interactome_file = 'interactom_nonredundant'
    re_output_file = 'filtered_file_v9.tsv'
    redundancy_remove.main(input_file, interactome_file ,re_output_file) 



################################################
# Autoencoder TRAINing
# with multiple inputs attributes

#Train atr
SAVE=True
batch_size =  32
pd_results = []
LOOP = True
processed_sample_values = [2000]
size_values = [7]
latent_dim_values = [2]    
epochs = [10]
ranges = 3

if A_RUN == True:
        for processed_sample in processed_sample_values:
                for latent_dim in latent_dim_values:
                    for size in size_values:
                        for epoch in epochs:
                            for i in range(ranges):
                                model_name = 'dina_model_sample_%d_dim_%d_size_%d_epochs_%d_index_%d.keras' % (processed_sample, latent_dim, size, epoch , i )
                                print(model_name)
                                print(f"Running with latent_dim={latent_dim}, epoch={epoch} , processed_sample_num={processed_sample} , matrix_size = {size} x {size}")
                                try:    
                                    encoded_vector, collected_data = data_autoencoding.main(latent_dim, model_name, processed_sample , size , SAVE , epoch, batch_size, auto_dir, pickle_ca_path, pickle_mean_path )
                                    #creating pandas dataframe for heatmap
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
                                    print(f"An IndexError occurred: {e}. Skipping to the next iteration.")
                                    continue
df = pd.DataFrame(pd_results)
tsv = 'autoencoder_train_sample_%d_dim_%d_size_%d_epochs_%d.tsv' % (processed_sample, latent_dim, size, epoch)
df.to_csv( tsv , sep='\t', index=False) 

################################################
# Autoencoder TESTing

#Test atr

pickle_path = '' # pickel_ca_path
if A_TEST == True:
    pickle_data = [os.path.join(pickle_dir, file) for file in os.listdir(pickle_dir) if file.endswith('.pickle')]
    for p in pickle_data:
        encode_test_model_name = ''  #'dina_model_sample_%d_dim_%d_size_%d_epochs_%d_index_%d.keras' % (processed_sample, latent_dim, size, epoch , i )
        print(encode_test_model_name)
        try:    
            encoded_vector_test = run_autoencoder.main(encode_test_model_name, sub_size, pickle_path)
        except IndexError as e:
            print(f"An IndexError occurred: {e}. Skipping to the next iteration.")

################################################
# Vector append

input_tsv_for_vec = '' #
vectorized_tsv_name = 'bert_test_with_vector_005'
adding_vector.main(pickle_dir, sub_size, input_tsv_for_vec , encode_test_model_name, vectorized_tsv_name )


################################################
#BERT Training

output_stat_tsv = 'BERT_training_stats_withvectors_004.tsv'
train_input_tsv = 'concatenated_file_v1.tsv'  #vectorized_tsv_name

if B_RUN == True :
    info = 'rep01020304'
    combined_fields = ["Residues","DSSP Structure","DSSP Index", "Vector"]
    try:
        full_df = pd.read_csv(output_stat_tsv, sep='\t')
    except FileNotFoundError:
        full_df = pd.DataFrame()

    for i in range(1, len(combined_fields) + 1):
        bert_model_name = 'bert_attrnum_%d_no_%s' % (i, info )
        print(f'Running model name : {bert_model_name}')
        fields_to_combine = combined_fields[:i]
        df_stats = combine_input_bert.main(bert_model_name, train_input_tsv, fields_to_combine)
        #tsv_output_path = model_name.replace('.pth', '_training_stats.tsv')
        #df_stats.to_csv('BERT_training_stats.tsv', sep='\t')
        df_stats.loc[0, 'model_name'] = bert_model_name
        full_df = pd.concat([full_df, df_stats], ignore_index=True)
    full_df.to_csv(output_stat_tsv, sep='\t', index=False)

################################################
#BERT Testing

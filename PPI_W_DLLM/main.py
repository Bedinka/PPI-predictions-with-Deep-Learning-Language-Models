import feature_extraction
import data_autoencoding
import combine_input_bert
import redundancy_remove

import pandas as pd
import os
F_RUN = True 
B_RUN = False
S_REMOVE = False
A_RUN = False
# Running feature extraction
# Creates the classes, the training .tsv , 

sample = 50
sub_size = 7
work_dir_f = "/home/dina/Documents/PPI_WDLLM/workdir"
tsv_f = "bert_small_test.tsv"
# Running Feature extraction
if F_RUN == True:
    feature_extraction.main(sample, sub_size, tsv_f)


# Removing non representative sequences

if S_REMOVE == True :
    input_file = tsv_f
    interactome_file = 'interactom_nonredundant'
    output_file = 'filtered_file_v9.tsv'
    redundancy_remove.main(input_file, interactome_file ,output_file)

# Autoencoder 

work_dir = "/home/dina/Documents/PPI_WDLLM"
processed_sample_values = [2000]
size_values = [7]
latent_dim_values = [2]    
epochs = [10]
SAVE=True
ranges = 3
TEST = False 
batch_size =  32
unseen_data_path =''
pd_results = []
if A_RUN == True:
    if TEST == False:
        for processed_sample in processed_sample_values:
                for latent_dim in latent_dim_values:
                    for size in size_values:
                        for epoch in epochs:
                            for i in range(ranges):
                                model_name = 'dina_model_sample_%d_dim_%d_size_%d_epochs_%d_index_%d.keras' % (processed_sample, latent_dim, size, epoch , i )
                                print(model_name)
                                print(f"Running with latent_dim={latent_dim}, epoch={epoch} , processed_sample_num={processed_sample} , matrix_size = {size} x {size}")
                                try:    
                                    encoded_vector, collected_data = data_autoencoding.main(latent_dim, model_name, processed_sample , size , SAVE , epoch , TEST, unseen_data_path, batch_size )
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
    else:
        pickle_data = [os.path.join(f"{work_dir}/Matrices_CA", file) for file in os.listdir(f"{work_dir}/Matrices_CA") if file.endswith('.pickle')]
        for p in pickle_data:
            unseen_data_path =  p
            processed_sample = 200
            latent_dim = 2
            size = 7
            epoch = 10
            i = 0
            model_name = 'dina_model_sample_%d_dim_%d_size_%d_epochs_%d_index_%d.keras' % (processed_sample, latent_dim, size, epoch , i )
            print(model_name)
            print(f"Running with latent_dim={latent_dim}, epoch={epoch} , processed_sample_num={processed_sample} , matrix_size = {size} x {size}")
            try:    
                encoded_vector, collected_data = data_autoencoding.main(latent_dim, model_name, processed_sample , size , SAVE , epoch , TEST, unseen_data_path, batch_size )
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

    df = pd.DataFrame(pd_results)
    tsv = 'train_sample_%d_dim_%d_size_%d_epochs_%d_index_%d.tsv' % (processed_sample, latent_dim, size, epoch , i )
    df.to_csv( tsv , sep='\t', index=False) 

# Vector append




#BERT 

out_csv = 'BERT_training_stats_withvectors_004.tsv'
tsv_path = 'concatenated_file_v1.tsv'
if B_RUN == True :
    info = 'rep01020304'
    combined_fields = ["Residues","DSSP Structure","DSSP Index", "Vector"]
    try:
        full_df = pd.read_csv(out_csv, sep='\t')
    except FileNotFoundError:
        full_df = pd.DataFrame()

    for i in range(1, len(combined_fields) + 1):
        bert_model_name = 'bert_attrnum_%d_no_%s' % (i, info )
        print(f'Running model name : {bert_model_name}')
        fields_to_combine = combined_fields[:i]
        df_stats = combine_input_bert.main(bert_model_name, tsv_path, fields_to_combine)
        #tsv_output_path = model_name.replace('.pth', '_training_stats.tsv')
        #df_stats.to_csv('BERT_training_stats.tsv', sep='\t')
        df_stats.loc[0, 'model_name'] = bert_model_name
        full_df = pd.concat([full_df, df_stats], ignore_index=True)
    full_df.to_csv(out_csv, sep='\t', index=False)


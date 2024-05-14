import feature_extraction
import data_autoencoding
import combine_input_bert

import pandas as pd

# Running feature extraction
# Creates the classes, the training .tsv , 
F_RUN = False
sample = 2000
sub_size = 7
work_dir_f = "/home/dina/Documents/PPI_WDLLM/workdir"
tsv_path = "bert_train_8.tsv"
# Running Feature extraction
if F_RUN == True:
    feature_extraction.main(sample, sub_size, tsv_path)

# Autoencoder 
A_RUN = False
work_dir = "/home/dina/Documents/PPI_WDLLM"
processed_sample_values = [1000]
size_values = [7]
latent_dim_values = [2]    
epochs = [10]
SAVE=True
ranges = 3
pd_results = []
for processed_sample in processed_sample_values:
        for latent_dim in latent_dim_values:
            for size in size_values:
                
                for epoch in epochs:
                    for i in range(ranges):

                        # Running Autoencoder
                        if A_RUN == True:
                            model_name = 'dina_model_sample_%d_dim_%d_size_%d_epochs_%d_index_%d.keras' % (processed_sample, latent_dim, size, epoch , i )
                            print(model_name)
                            #fout.write(model_path+"\n")
                            print(f"Running with latent_dim={latent_dim}, epoch={epoch} , processed_sample_num={processed_sample} , matrix_size = {size} x {size}")
                            try:    
                                collected_data = data_autoencoding.main(latent_dim, model_name, processed_sample , size , SAVE , epoch )
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
tsv = 'train_sample_%d_dim_%d_size_%d_epochs_%d_index_%d.tsv' % (processed_sample, latent_dim, size, epoch , i )
df.to_csv( tsv , sep='\t', index=False) 

#BERT 
B_RUN = True
if B_RUN == True :
    tsv_num = 8
    tsv_path = 'bert_train_%d.tsv' % (tsv_num)
    combined_fields = ["Protein ID","Residues","DSSP Structure","DSSP Index"]
    try:
        full_df = pd.read_csv('BERT_training_stats_v1.tsv', sep='\t')
    except FileNotFoundError:
        full_df = pd.DataFrame()

    for i in range(1, len(combined_fields) + 1):
        bert_model_name = 'bert_model_attrnum_%d_tsvnum_%d.pth' % (i, tsv_num )
        print(f'Running model name : {bert_model_name}')
        fields_to_combine = combined_fields[:i]
        df_stats = combine_input_bert.main(bert_model_name, tsv_path, fields_to_combine)
        # Save training stats to TSV
        #tsv_output_path = model_name.replace('.pth', '_training_stats.tsv')
        #df_stats.to_csv('BERT_training_stats.tsv', sep='\t')
        df_stats.loc[0, 'model_name'] = bert_model_name
        # Append the new statistics to the full dataframe
        full_df = pd.concat([full_df, df_stats], ignore_index=True)

    # Save the full DataFrame back to the TSV file
    full_df.to_csv('BERT_training_stats_v1.tsv', sep='\t', index=False)


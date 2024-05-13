import feature_extraction
import data_autoencoding
import combine_input_bert

import pandas as pd

# Running feature extraction
# Creates the classes, the training .tsv , 
F_RUN = False  
processed_sample_values = [100]
size_values = [7]
work_dir_f = "/home/dina/Documents/PPI_WDLLM/workdir"

# Autoencoder 
A_RUN = False
work_dir = "/home/dina/Documents/PPI_WDLLM"
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
                        # Running Feature extraction
                        if F_RUN == True:
                             feature_extraction.main(processed_sample, size)
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
model_name = 'bert_model_v7.pth'
tsv_path = 'bert_train_6.tsv'
combine_input_bert.main(model_name, tsv_path) 

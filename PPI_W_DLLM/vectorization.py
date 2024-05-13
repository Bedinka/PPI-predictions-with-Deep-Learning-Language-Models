import data_autoencoding    
import pandas as pd

def main():
    
    processed_sample_values = [500]  
    latent_dim_values = [2]  
    size_values = [7]  
    epochs = [10]
    results = {}
    pd_results = []
    SAVE=True
    for processed_sample in processed_sample_values:
        for latent_dim in latent_dim_values:
            for size in size_values:
                for epoch in epochs:
                    for i in range(3):
                        model_name = 'dina_model_sample_%d_dim_%d_size_%d_epochs_%d_index_%d.keras' % (processed_sample, latent_dim, size, epoch , i )
                        print(model_name)
                        #fout.write(model_path+"\n")
                        print(f"Running with latent_dim={latent_dim}, epoch={epoch} , processed_sample_num={processed_sample} , matrix_size = {size} x {size}")
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
    df = pd.DataFrame(pd_results)
    tsv = 'train_sample_%d_dim_%d_size_%d_epochs_%d_index_%d.tsv' % (processed_sample, latent_dim, size, epoch , i )
    df.to_csv( tsv , sep='\t', index=False) 



if __name__ == "__main__":
    main()
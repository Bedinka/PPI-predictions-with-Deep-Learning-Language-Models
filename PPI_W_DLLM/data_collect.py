import data_autoencoding
import pandas as pd

def main():
    
    processed_sample_values = [5, 10]  
    size_values = [7 , 8, 9 ]  
    latent_dim_values = [2, 4, 8]  
    epochs = [ 10, 20]
    results = {}
    pd_results = []
    
    SAVE=True
    for processed_sample in processed_sample_values:
        for latent_dim in latent_dim_values:
            for size in size_values:
                for epoch in epochs:
                    fout = open("summary_detailed_3.txt", "a")
                    for i in range(2):
                        model_name = 'dina_model_sample_%d_dim_%d_size_%d_epochs_%d_index_%d.keras' % (processed_sample, latent_dim, size, epoch , i )
                        print(model_name)
                        #fout.write(model_path+"\n")
                        print(f"Running with latent_dim={latent_dim}, epoch={epoch} , processed_sample_num={processed_sample} , matrix_size = {size} x {size}")
                        collected_data = data_autoencoding.main(latent_dim, model_name, processed_sample , size , SAVE , epoch )
                        #creating txt 
                        fout.write(str(collected_data)+"\n"+"\n")
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
    fout.close()
    df = pd.DataFrame(results)
    df.to_csv("results.csv", index=False)

    return df


if __name__ == "__main__":
    main()
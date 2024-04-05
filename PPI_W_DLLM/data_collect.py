import data_autoencoding

def main():
    
    processed_sample_values = [5, 10, 15 , 20 ]  
    size_values = [7, 8, 9 , 10 , 11 ,12]  
    latent_dim_values = [2, 4, 7]  
    model_path_values = ['dina_model.keras']  
    epochs = 10
    results = {}
    for latent_dim in latent_dim_values:
        for model_path in model_path_values:
            for processed_sample in processed_sample_values:
                for size in size_values:
                    print(f"Running with latent_dim={latent_dim}, model_path={model_path} , processed_sample_num={processed_sample}, matrix_size = {size} x {size}")
                    losses = data_autoencoding.main(latent_dim, model_path, processed_sample , size ,epochs )
                    results[(latent_dim, model_path)] = losses
    

    print(results)

if __name__ == "__main__":
    main()
import feature_extraction
import data_autoencoding

def main():
    
    processed_sample_values = [5, 10, 15]  
    size_values = [7, 10, 15]  
    latent_dim_values = [2, 5, 10]  
    model_path_values = ['model1.keras', 'model2.keras', 'model3.keras']  
    epochs = 10

    for processed_sample in processed_sample_values:
        for size in size_values:
            print(f"Running with processed_sample={processed_sample} and size={size}")
            interacting_proteins = feature_extraction.main(processed_sample, size)
            # Process the obtained interacting_proteins as needed
    results = {}
    for latent_dim in latent_dim_values:
        for model_path in model_path_values:
            print(f"Running with latent_dim={latent_dim} and model_path={model_path}")
            losses = data_autoencoder.main(latent_dim, model_path, epochs)
            results[(latent_dim, model_path)] = losses

    print(results)

if __name__ == "__main__":
    main()
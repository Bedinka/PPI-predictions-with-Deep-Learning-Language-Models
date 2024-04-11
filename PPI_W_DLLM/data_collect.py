import data_autoencoding

def main():
    
    processed_sample_values = [5, 10, 15 , 20 ]  
    size_values = [7, 8, 9 , 10 , 11 ,12]  
    latent_dim_values = [2, 3, 4, 5, 6, 7]  
    epochs = [ 10, 20, 40]
    results = {}
    
    SAVE=True
    for latent_dim in latent_dim_values:
        for processed_sample in processed_sample_values:
            for size in size_values:
                for epoch in epochs:
                    fout = open("summary.txt", "a")
                    for i in range(10):
                        model_path = 'dina_model_sample_%d_dim_%d_size_%d_index_%d.keras' % (processed_sample, latent_dim, size, i)
                        print(model_path)
                        fout.write(model_path+"\n")
                        print(f"Running with latent_dim={latent_dim}, epoch={epoch} , processed_sample_num={processed_sample}, matrix_size = {size} x {size}")
                        collected_data = data_autoencoding.main(latent_dim, model_path, processed_sample , size , SAVE , epoch )
                        print(collected_data)
                        fout.write(str(collected_data)+"\n"+"\n")  
    fout.close()


    print(results)

if __name__ == "__main__":
    main()
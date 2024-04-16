import data_autoencoding

def main():
    
    processed_sample_values = [5, 15]  
    size_values = [7 , 8, 9 ]  
    latent_dim_values = [2, 4, 8]  
    epochs = [ 10, 20]
    results = {}
    
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
                        #print(collected_data)
                        fout.write(str(collected_data)+"\n"+"\n")  
    fout.close()


if __name__ == "__main__":
    main()
import numpy as np
import YangLabIntern.PPI_W_DLLM.data_autoencoding as data_autoencoding
from transformers import BertConfig, BertTokenizer, BertModel
config = BertConfig.from_pretrained('bert-base-cased')
tokenizer = BertTokenizer.from_pretrained('bert-base-cased')
model = BertModel.from_pretrained('bert-base-cased', config=config)

def main():
    latent_dim = 4
    processed_sample=1
    size=7
    SAVE = True
    epoch =10
    i=1
    model_name= 'dina_model_sample_%d_dim_%d_size_%d_epochs_%d_index_%d.keras' % (processed_sample, latent_dim, size, epoch , i )
    processed_sample=5

    # Example numerical matrix input
    numerical_matrix = np.array(data_autoencoding.main(latent_dim, model_name, processed_sample, size, SAVE, epoch))
    print(numerical_matrix.shape)
    # Convert numerical matrix to string representation
    text_data = [[str(num) for num in row] for row in numerical_matrix]

    # Flatten the list of lists into a single list
    flat_text_data = [item for sublist in text_data for item in sublist]

    tokenizer = BertTokenizer.from_pretrained('bert-base-uncased')

    # Tokenize the data
    tokenized_data = [tokenizer.tokenize(token) for token in flat_text_data]

    print("Data tokenized")

    return tokenized_data , numerical_matrix

if __name__ == "__main__":
    main()
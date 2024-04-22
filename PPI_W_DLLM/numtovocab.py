import numpy as np
import data_autoencoding
import torch
from transformers import BertTokenizer, BertForSequenceClassification, AdamW
from torch.utils.data import DataLoader, Dataset


def fine_tune_bert(autoencoder_output, labels):
    # Tokenize autoencoder output (if needed)
    # Here we just use it as is, but if it's not textual, you'll need to come up with tokenization logic
    tokenized_data = autoencoder_output

    # Convert data and labels to PyTorch tensors
    tokenized_data = torch.tensor(tokenized_data, dtype=torch.float32)
    labels = torch.tensor(labels, dtype=torch.long)

    # Define a custom dataset
    class CustomDataset(Dataset):
        def __init__(self, data, labels):
            self.data = data
            self.labels = labels

        def __len__(self):
            return len(self.data)

        def __getitem__(self, idx):
            input_data = self.data[idx % len(self.data)]  # Modulo operation to handle out-of-bounds index
            label = self.labels[idx % len(self.labels)]   # Modulo operation to handle out-of-bounds index
            return {'input': input_data, 'labels': label}

    # Create dataset and dataloader
    dataset = CustomDataset(tokenized_data, labels)
    dataloader = DataLoader(dataset, batch_size=4, shuffle=True)

    # Load pre-trained BERT model and tokenizer
    tokenizer = BertTokenizer.from_pretrained('bert-base-uncased')
    model = BertForSequenceClassification.from_pretrained('bert-base-uncased', num_labels=2)

    # Fine-tuning BERT
    optimizer = AdamW(model.parameters(), lr=1e-5)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)

    model.train()
    for epoch in range(3):  # Adjust number of epochs as needed
        for batch in dataloader:
            optimizer.zero_grad()
            input_ids = batch['input'].to(device).long()
            labels = batch['labels'].to(device)
            outputs = model(input_ids, labels=labels)
            loss = outputs.loss
            loss.backward()
            optimizer.step()

# Example usage
def main():
    latent_dim = 2
    processed_sample = 5
    size = 7
    SAVE = True
    epoch = 10
    i = 1
    model_name = 'dina_model_sample_%d_dim_%d_size_%d_epochs_%d_index_%d.keras' % (processed_sample, latent_dim, size, epoch, i)
    processed_sample = 5

    # Example numerical matrix input
    numerical_matrix = np.array(data_autoencoding.main(latent_dim, model_name, processed_sample, size, SAVE, epoch))
    print(numerical_matrix.shape)


    # Example labels
    labels = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]  

    # Fine-tune BERT
    fine_tune_bert(numerical_matrix, labels)

if __name__ == "__main__":
    main()

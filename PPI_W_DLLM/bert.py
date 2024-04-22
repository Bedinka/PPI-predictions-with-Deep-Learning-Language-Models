import torch
from transformers import BertTokenizer, BertForSequenceClassification
from torch.utils.data import DataLoader, Dataset
from transformers import AdamW
import numtovocab

# Sample data
autoencoder_output = numtovocab.main()  # Example autoencoder output

# Tokenize autoencoder output (if needed)
# Here we just use it as is, but if it's not textual, you'll need to come up with tokenization logic
tokenized_data = autoencoder_output

# Prepare labels (if any)
labels = [0, 1, 1, ...]  # Example labels

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
        return {'input': self.data[idx], 'labels': self.labels[idx]}

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
        input_ids = batch['input'].to(device)
        labels = batch['labels'].to(device)
        outputs = model(input_ids, labels=labels)
        loss = outputs.loss
        loss.backward()
        optimizer.step()

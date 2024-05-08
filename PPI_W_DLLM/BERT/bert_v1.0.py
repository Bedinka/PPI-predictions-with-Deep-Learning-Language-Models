# My BERT try fith flattened str data mfrom numerical input 

import torch
import numpy as np 
from transformers import BertTokenizer, BertForSequenceClassification
from torch.utils.data import DataLoader, Dataset
from transformers import AdamW
import numtovocab

tokenized , numerical_matrix = numtovocab.main() 

labels = np.random.randint(0, 2, size=100)
tokenized_data = torch.tensor(numerical_matrix, dtype=torch.long)
labels = torch.tensor(labels, dtype=torch.long)

class CustomDataset(Dataset):
    def __init__(self, numerical_matrix, labels):
        self.numerical_matrix = numerical_matrix
        self.labels = labels

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        return {'input': self.numerical_matrix[idx], 'labels': self.labels[idx]}


dataset = CustomDataset(tokenized_data, labels)
dataloader = DataLoader(dataset, batch_size=4, shuffle=True)

tokenizer = BertTokenizer.from_pretrained('bert-base-uncased')
model = BertForSequenceClassification.from_pretrained('bert-base-uncased', num_labels=2)

optimizer = AdamW(model.parameters(), lr=1e-5)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

model.train()

for epoch in range(3):  
    for batch in dataloader:
        optimizer.zero_grad()
        input_ids = batch['input'].to(device)
        labels = batch['labels'].to(device)
        outputs = model(input_ids, labels=labels)
        loss = outputs.loss
        loss.backward()
        optimizer.step()

torch.save(model.state_dict(), 'bert_model_state_dict.pth')
print("Model trained and saved")
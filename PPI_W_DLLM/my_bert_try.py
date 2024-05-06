import torch
from transformers import BertForSequenceClassification, BertTokenizer
import numtovocab
# Sample data (you can replace this with your actual data)
tokenized , numerical_matrix = numtovocab.main()
input_text = str(numerical_matrix)
labels = [0, 1]

# Tokenize input text
tokenizer = BertTokenizer.from_pretrained('bert-base-uncased')
tokenized_input = tokenizer(input_text, padding=True, truncation=True, return_tensors='pt')

# Define and train BERT model (you can replace this with your actual training code)
model = BertForSequenceClassification.from_pretrained('bert-base-uncased', num_labels=2)
optimizer = torch.optim.AdamW(model.parameters(), lr=1e-5)

for epoch in range(3):  # Adjust number of epochs as needed
    for batch in range(len(input_text)):
        optimizer.zero_grad()
        output = model(**tokenized_input, labels=torch.tensor([labels[batch]]))
        loss = output.loss
        loss.backward()
        optimizer.step()

model_save_path = 'bert_model_state_dict.pth'
torch.save({
    'state_dict': model.state_dict(),
    'tokenizer': tokenizer,
    'config': model.config
}, model_save_path)

loaded_model_dict = torch.load(model_save_path)
loaded_model = BertForSequenceClassification.from_pretrained('bert-base-uncased', num_labels=2)
loaded_model.load_state_dict(loaded_model_dict['state_dict'])


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
loaded_model.to(device)

# Example inference with the loaded model
input_text = ["New input text for inference"]
tokenized_input = tokenizer(input_text, padding=True, truncation=True, return_tensors='pt').to(device)
with torch.no_grad():
    outputs = loaded_model(**tokenized_input)
    predictions = torch.argmax(outputs.logits, dim=1)

print("Predictions:", predictions)

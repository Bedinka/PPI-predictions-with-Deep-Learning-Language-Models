import torch
import pandas as pd
from transformers import BertTokenizer, BertForSequenceClassification

tsv_path = 'ftest_with_vector_ca.tsv'
data_df = pd.read_csv( tsv_path , sep='\t')
data_df.head()
sen_w_feats = []
labels = []
data_df = data_df.fillna("")

# Define a random sequence
combined_fields = ["Residues","DSSP Structure","DSSP Index", "Vector"]

for index, row in data_df.iterrows():

        fields = [str(row[field]) for field in combined_fields]
        combined = '[SEP]'.join(fields)
        sen_w_feats.append(combined)
        labels.append(row["Interact"])

tokenizer = BertTokenizer.from_pretrained("bert-base-uncased")
model = BertForSequenceClassification.from_pretrained("bert")

inputs = tokenizer(combined, return_tensors="pt")

with torch.no_grad():
    outputs = model(**inputs)

predicted_label = torch.argmax(outputs.logits).item()

print("Predicted label:", predicted_label)

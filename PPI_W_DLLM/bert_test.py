import torch
import pandas as pd
from transformers import BertTokenizer, BertForSequenceClassification
from sklearn.metrics import roc_curve, auc, accuracy_score
import matplotlib.pyplot as plt

# Load data
tsv_path = 'bert_test_with_vector_004.tsv'
data_df = pd.read_csv(tsv_path, sep='\t')
data_df = data_df.fillna("")

sen_w_feats = []
labels = []
combined_fields = ["Residues", "DSSP Structure", "DSSP Index", "Vector"]

for index, row in data_df.iterrows():
    fields = [str(row[field]) for field in combined_fields]
    combined = '[SEP]'.join(fields)
    sen_w_feats.append(combined)
    labels.append(row["Interact"])  

# Load tokenizer and model
tokenizer = BertTokenizer.from_pretrained("bert-base-uncased")
model_path = 'bert_attrnum_4_no_rep01020304'  
model = BertForSequenceClassification.from_pretrained(model_path)

def process_sequence(sequence):
    encoded_dict = tokenizer.encode_plus(
        sequence,
        add_special_tokens=True,
        max_length=512,
        truncation=True,
        padding='max_length',
        return_attention_mask=True,
        return_tensors='pt',
    )
    return encoded_dict

# Move model to appropriate device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

# Make predictions
predicted_labels = []
logit_list = []
for seq in sen_w_feats:
    encoded_dict = process_sequence(seq)
    input_ids = encoded_dict['input_ids'].to(device)
    attention_mask = encoded_dict['attention_mask'].to(device)

    with torch.no_grad():
        outputs = model(input_ids, token_type_ids=None, attention_mask=attention_mask)
        logits = outputs.logits
    predicted_label = torch.argmax(logits).item()
    predicted_labels.append(predicted_label)
    logit_list.append(logits.cpu().numpy())

# Calculate accuracy
accuracy = accuracy_score(labels, predicted_labels)
print(f"Accuracy: {accuracy:.4f}")

# Calculate misclassifications
misclassified_as_0 = sum(1 for predicted, actual in zip(predicted_labels, labels) if predicted == 0 and actual == 1)
misclassified_as_1 = sum(1 for predicted, actual in zip(predicted_labels, labels) if predicted == 1 and actual == 0)
total_1s = sum(1 for label in labels if label == 1)
total_0s = sum(1 for label in labels if label == 0)

percent_1_labeled_as_0 = (misclassified_as_0 / total_1s) * 100 if total_1s != 0 else 0
percent_0_labeled_as_1 = (misclassified_as_1 / total_0s) * 100 if total_0s != 0 else 0

print(f"Percentage of 1s labeled as 0: {percent_1_labeled_as_0:.2f}%")
print(f"Percentage of 0s labeled as 1: {percent_0_labeled_as_1:.2f}%")

# Calculate ROC Curve and AUC
logit_list = [logit[0] for logit in logit_list]  # Remove unnecessary dimensions
logit_scores = [logit[1] for logit in logit_list]  # Get the logit scores for the positive class

fpr, tpr, _ = roc_curve(labels, logit_scores)
roc_auc = auc(fpr, tpr)

print(f"AUC: {roc_auc:.4f}")

# Plot ROC Curve
plt.figure()
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc="lower right")
plt.show()

# Print individual predictions
'''for predicted_label, actual_label in zip(predicted_labels, labels):
    print(f"Predicted label: {predicted_label}, Actual label: {actual_label}")
'''
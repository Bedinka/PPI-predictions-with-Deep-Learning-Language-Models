import torch
import pandas as pd
import os
import datetime
import logging
import warnings
from transformers import BertTokenizer, BertForSequenceClassification
from sklearn.metrics import roc_curve, auc, accuracy_score
import matplotlib.pyplot as plt
import setup_run

import os
import sys
import datetime
import logging
import itertools

class LoggerWriter:
    def __init__(self, level):
        self.level = level

    def write(self, message):
        if message.strip():
            self.level(message)

    def flush(self):
        pass

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

def load_vectors_for_bert(input_tsv, vector_pickle_dir):
    df = pd.read_csv(input_tsv, sep='\t')
    vectors = []
    import pickle
    for index, row in df.iterrows():
        vector_pickle_path = os.path.join(vector_pickle_dir, f'vector_{row["Protein ID"]}.pickle')
        if os.path.exists(vector_pickle_path):
            with open(vector_pickle_path, 'rb') as f:
                vector = pickle.load(f)
            vectors.append(vector)
        else:
            vectors.append(None)  # or handle missing vectors appropriately
    
    df['Vector'] = vectors
    return df

set = False 
path = './2024-06-13_15-03-02'
run_directory = setup_run.setup_run_directory(set, path)

# Logger 
time = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')


# Load data
vector_pickle_dir = '/home/dina/Documents/PPI_WDLLM/Matrices_CA/Vector_pickle_test'
tsv_path = '/home/dina/Documents/PPI_WDLLM/2024-06-13_13-00-08/2024-06-13_13-00-08_TEST_dataset.tsv'
#data_df = pd.read_csv(tsv_path, sep='\t')
data_df = load_vectors_for_bert(tsv_path, vector_pickle_dir)
data_df = data_df.fillna("")

'''# Model choice
combined_fields = ["Residues", "DSSP Structure", "DSSP Index", "RSA", "Vector"]
logging.info('Testing BERT with %s , in path %s and with pickle load' , tsv_path, path)

for i in range(1, len(combined_fields) + 1):
    fields_to_combine = combined_fields[:i]
    sen_w_feats = [] 
    labels = []
    bert_dir = '/home/dina/Documents/PPI_WDLLM/2024-06-11_12-36-28/BERT/'
    model_name = 'bert_attrnum_%d_no_1300_s_wDataloader_interactom'  %(i)
    logging.info( "#" * 50 )
    logging.info( "# BERT TESTING \n" "%s" , model_name )
    
    logging.info(" With combined fields : %s " , fields_to_combine)
    for index, row in data_df.iterrows():
        fields = [str(row[field]) for field in fields_to_combine]
        combined = '[SEP]'.join(fields)
        sen_w_feats.append(combined)
        labels.append(row["Interact"])  '''
# Model choice
input_fields = ["Residues", "DSSP Structure", "DSSP Index", "RSA", "Vector"]
logging.info('Testing BERT with %s, in path %s and with pickle load', tsv_path, vector_pickle_dir)
# Generate all combinations of input fields
field_combinations = []
for r in range(1, len(input_fields) + 1):
    field_combinations.extend(itertools.combinations(input_fields, r))

for fields_to_combine in field_combinations:
    fields_to_combine = list(fields_to_combine)
    sen_w_feats = []
    labels = []
    bert_dir = '/home/dina/Documents/PPI_WDLLM/2024-06-13_15-03-02/BERT/'
    model_name = 'bert_combined_fields_%s_%s' % ('1500s_allcombinations', '_'.join(fields_to_combine))
    logging.info("#" * 50)
    logging.info("# BERT TESTING \n%s", model_name)
    logging.info(" With combined fields: %s", fields_to_combine)
    
    for index, row in data_df.iterrows():
        fields = [str(row[field]) for field in fields_to_combine]
        combined = '[SEP]'.join(fields)
        sen_w_feats.append(combined)
        labels.append(row["Interact"])

    # Load tokenizer and model
    tokenizer = BertTokenizer.from_pretrained("bert-base-uncased")
    directory = os.path.join(bert_dir, model_name)
    print(model_name)
    model = BertForSequenceClassification.from_pretrained(directory)
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

    # Calculate true positives, true negatives, false positives, and false negatives
    tp = sum(1 for predicted, actual in zip(predicted_labels, labels) if predicted == 1 and actual == 1)
    tn = sum(1 for predicted, actual in zip(predicted_labels, labels) if predicted == 0 and actual == 0)
    fp = sum(1 for predicted, actual in zip(predicted_labels, labels) if predicted == 1 and actual == 0)
    fn = sum(1 for predicted, actual in zip(predicted_labels, labels) if predicted == 0 and actual == 1)

    print(f"True Positives (TP): {tp}")
    print(f"True Negatives (TN): {tn}")
    print(f"False Positives (FP): {fp}")
    print(f"False Negatives (FN): {fn}")

    # Calculate sensitivity and specificity
    sensitivity = tp / (tp + fn) if (tp + fn) != 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) != 0 else 0

    print(f"Sensitivity: {sensitivity:.4f}")
    print(f"Specificity: {specificity:.4f}")

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
    plt.savefig(model_name)

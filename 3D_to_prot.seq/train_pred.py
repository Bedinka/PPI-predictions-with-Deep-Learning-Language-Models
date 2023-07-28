import csv
from transformers import BertTokenizerFast,BertForMaskedLM, BertTokenizer, pipeline
import json

 # Read JSON config file
with open("./config/config.json", "r") as file:
    config = json.load(file)

tokenizer_path = config["tokenizer_path"]
tokenizer = BertTokenizerFast.from_pretrained(tokenizer_path)
model = BertForMaskedLM.from_pretrained("/home/oem/Desktop/ALEX/Project/proteinbert/model")
unmasker = pipeline('fill-mask', model=model, tokenizer=tokenizer, top_k=1)

tokenizer.add_tokens("J")
model.resize_token_embeddings(len(tokenizer))

# Load sequences from FASTA file
sequences = {}
with open("/home/oem/Desktop/ALEX/Project/data/structure10_250.fa") as f:
    sequence_id = ""
    sequence = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if sequence_id:
                sequences[sequence_id] = sequence
                print(sequence_id)
            sequence_id = line[1:]
            sequence = ""
        else:
            sequence += line
    sequences[sequence_id] = sequence
    print(sequence_id)

# Predict missing letters for each sequence and write to CSV file
MOTIF_LEN = 20
with open("score2.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["ID", "Sequence"])
    
    for sequence_id, sequence in sequences.items():
        #predicted_sequence= ""
        #predicted_sequence1 = ""
        predicted_sequence2 = ""
        #for i in range(MOTIF_LEN, MOTIF_LEN+20): #len(sequence)-MOTIF_LEN):
        for i in range(MOTIF_LEN, len(sequence)-MOTIF_LEN):
            #masked_sequence = " ".join(sequence[i-10:i]) + " [MASK] " + " ".join(sequence[i+1:i+11])
            #masked_sequence1 = " ".join(sequence[i-MOTIF_LEN:i]) + " [MASK] " + " ".join(sequence[i+1:i+1+MOTIF_LEN])
            masked_sequence2 = " ".join(sequence[:i]) + " [MASK] " + " ".join(sequence[i+1:])
            #masked_sequence = sequence[:i] + "[MASK]" + sequence[i+1:]
            #prediction1 = unmasker(masked_sequence1) #[0]["token_str"]
            #prediction= unmasker(masked_sequence)
            prediction2 = unmasker(masked_sequence2) #[0]["token_str"]
            #predicted_sequence += prediction[0]["token_str"]
            #predicted_sequence1 += prediction1[0]["token_str"]
            predicted_sequence2 += prediction2[0]["token_str"]
            #print(predicted_sequence)
            #predicted_sequence2= predicted_sequence2.lower()
            #predicted_sequence2= predicted_sequence2.replace('z', 'Z')
            print(predicted_sequence2)
            #writer.writerow([sequence_id, sequence[i], prediction1[0]["token_str"], prediction1[0]["score"], prediction2[0]["token_str"], prediction2[0]["score"]])   
        #writer.writerow([sequence_id, predicted_sequence1])
        writer.writerow([sequence_id, predicted_sequence2])
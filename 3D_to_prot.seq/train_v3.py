import sys
import random
from sklearn.model_selection import train_test_split
from tqdm import tqdm
from torch.utils.data import DataLoader
import torch
from torch.utils.data import Dataset
from transformers import (
    BertTokenizerFast,
    BertForMaskedLM,
    BertConfig,
    TrainingArguments,
    Trainer,
    logging,
    pipeline
)
import os
import json



# Set a random seed for reproducibility
random.seed(42)

original_stdout = sys.stdout
with open('output.txt', 'w') as f:
    sys.stdout = f

    # Read JSON config file
    with open("./config/config.json", "r") as file:
        config = json.load(file)

    # Initial parameters
    os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:1024"
    torch.backends.cuda.matmul.allow_tf32 = True

    # Directory settings
    INPUT_DIR = config["input_dir"]
    OUTPUT_DIR = config["output_dir"]
    OUTPUT_FILE = config["output_file"]


    # Read input and output fasta sequences from .txt files
    def read_fasta_file(file_path):
        sequences = []
        sequence_ids = []  # List to store sequence IDs
        with open(file_path, "r") as file:
            lines = file.readlines()
            sequence = []
            sequence_id = None  # Variable to store the current sequence ID
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if sequence:
                        sequences.append("".join(sequence).upper())
                        sequence_ids.append(sequence_id)  # Append the current sequence ID
                    sequence = []
                    sequence_id = line[1:]  # Extract the ID, excluding the '>'
                else:
                    sequence.append(line)
            if sequence:
                sequences.append("".join(sequence).upper())
                sequence_ids.append(sequence_id)  # Append the last sequence ID
        return sequences, sequence_ids  # Return both the sequences and their corresponding IDs

    input_sequences, input_sequence_ids = read_fasta_file(INPUT_DIR)
    output_sequences, output_sequence_ids = read_fasta_file(OUTPUT_DIR) 
    input_train_sequences, input_val_sequences, input_train_sequence_ids, input_val_sequence_ids = train_test_split(input_sequences, input_sequence_ids, test_size=0.1, random_state=42)
    output_train_sequences, output_val_sequences, output_train_sequence_ids, output_val_sequence_ids = train_test_split(output_sequences, output_sequence_ids, test_size=0.1, random_state=42)

     
    # Model settings
    tokenizer_path = config["tokenizer_path"]
    tokenizer = BertTokenizerFast.from_pretrained(tokenizer_path)

    model_path = config["model_path"]
    model = BertForMaskedLM.from_pretrained(model_path).to("cuda")
    tokenizer.add_tokens("J")
    model.resize_token_embeddings(len(tokenizer))
    
    # Train arguments
    logging.set_verbosity_error()
    training_args = TrainingArguments(
        output_dir=OUTPUT_FILE,
        **config["training_arguments"],
        num_train_epochs=config["num_epochs"],
        disable_tqdm=False
    )   

    # Define dataset class
    class CustomDataset(Dataset):
        def __init__(self, sequences, labels, tokenizer, max_length):
            self.sequences = sequences
            self.labels = labels
            self.tokenizer = tokenizer
            self.max_length = max_length
    

        def __len__(self):
            return len(self.sequences)

        def __getitem__(self, idx):
            sequence = self.sequences[idx]
            label = self.labels[idx]
            sequence = sequence[:self.max_length]
            label = label[:self.max_length] 

            # Convert sequence and label in a str
            sequence = " ".join(sequence)
            label = " ".join(label)

            #PAD de small sequence
            padding_length = self.max_length - len(sequence)
            sequence = sequence + " [PAD]" * padding_length
            label = label + " [PAD]" * padding_length

            encoded_inputs = self.tokenizer(
                sequence,
                padding='max_length',
                truncation=True,
                max_length=self.max_length,
                return_tensors="pt",
            )
            encoded_labels = self.tokenizer(
                label,
                padding='max_length',
                truncation=True,
                max_length=self.max_length,
                return_tensors="pt",
            )
            labels = [
                -100 if token == self.tokenizer.pad_token_id else encoded_labels["input_ids"][0][i]
                for i, token in enumerate(encoded_inputs["input_ids"][0])
            ]
            
            return {
                "input_ids": encoded_inputs["input_ids"][0],
                "attention_mask": encoded_inputs["attention_mask"][0],
                "labels": torch.tensor(labels),
            }
    
    # Create dataset from sequences
    max_length = 250
    train_dataset = CustomDataset(input_train_sequences, output_train_sequences, tokenizer, max_length)
    val_dataset = CustomDataset(input_val_sequences, output_val_sequences, tokenizer, max_length)

    # Create dataloader
    batch_size = config["batch_size"]
    input_dataloader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    

    # Define trainer
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_dataset,
        eval_dataset=val_dataset,
        data_collator=lambda data: {
            "input_ids": torch.stack([item["input_ids"] for item in data]),
            "attention_mask": torch.stack([item["attention_mask"] for item in data]),
            "labels": torch.stack([item["labels"] for item in data]),
        }
    )

    # Trainer training and model saving
    trainer.train()
    progress_bar = tqdm(input_dataloader, desc=f"Processing Epoch 10")
    model.save_pretrained("./model")

sys.stdout = original_stdout

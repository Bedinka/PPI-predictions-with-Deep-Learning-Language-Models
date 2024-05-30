import pandas as pd
import torch
import os
from torch.utils.data import TensorDataset, DataLoader, SequentialSampler

        
def main(model_name, tsv_path, combined_fields): 
    if os.path.isfile(model_name):
        print(f"Loading saved model from {model_name}...")
        from transformers import BertTokenizer, BertForSequenceClassification
        import numpy as np
        # Load the BERT tokenizer and the model
        tokenizer = BertTokenizer.from_pretrained('bert-base-uncased', do_lower_case=True)
        model = BertForSequenceClassification.from_pretrained("bert-base-uncased", num_labels=2)
        model.load_state_dict(torch.load(model_name))
        model.cuda()
        
        # Load data for evaluation
        data_df = pd.read_csv(tsv_path, sep='\t')
        data_df = data_df.fillna("")
        
        # Combine fields into strings
        sen_w_feats = []
        labels = []
        
        for index, row in data_df.iterrows():
            fields = [str(row[field]) for field in combined_fields]
            combined = '[SEP]'.join(fields)
            sen_w_feats.append(combined)
            labels.append(row["Interact"])
        
        # Tokenize and encode the sentences
        input_ids = []
        attention_masks = []
        
        for sent in sen_w_feats:
            encoded_dict = tokenizer.encode_plus(
                                sent,                     
                                add_special_tokens=True,
                                max_length=512,           
                                truncation=True,
                                padding='max_length',
                                return_attention_mask=True,
                                return_tensors='pt',
                        )
            input_ids.append(encoded_dict['input_ids'])
            attention_masks.append(encoded_dict['attention_mask'])
        
        input_ids = torch.cat(input_ids, dim=0)
        attention_masks = torch.cat(attention_masks, dim=0)
        labels = torch.tensor(labels)
        
        # Split the dataset into train, validation, and test sets
        
        train_size = int(0.8 * len(data_df))
        val_size = int(0.1 * len(data_df))
        test_size = len(data_df) - (train_size + val_size)
        
        indices = np.arange(0, len(data_df))
        np.random.shuffle(indices)
        
        train_idx = indices[0:train_size]
        val_idx = indices[train_size:(train_size + val_size)]
        test_idx = indices[(train_size + val_size):]
        
        train_dataset = TensorDataset(input_ids[train_idx], attention_masks[train_idx], labels[train_idx])
        val_dataset = TensorDataset(input_ids[val_idx], attention_masks[val_idx], labels[val_idx])
        test_dataset = TensorDataset(input_ids[test_idx], attention_masks[test_idx], labels[test_idx])
        
        train_dataloader = DataLoader(train_dataset, sampler=SequentialSampler(train_dataset), batch_size=8)
        validation_dataloader = DataLoader(val_dataset, sampler=SequentialSampler(val_dataset), batch_size=8)
        prediction_dataloader = DataLoader(test_dataset, sampler=SequentialSampler(test_dataset), batch_size=8)
        
        # Evaluate the model
        def flat_accuracy(preds, labels):
            pred_flat = np.argmax(preds, axis=1).flatten()
            labels_flat = labels.flatten()
            return np.sum(pred_flat == labels_flat) / len(labels_flat)
        
        import time
        import datetime
        from sklearn.metrics import f1_score
        
        def format_time(elapsed):
            elapsed_rounded = int(round((elapsed)))
            return str(datetime.timedelta(seconds=elapsed_rounded))
        
        print("Predicting labels for test sentences...")
        
        model.eval()
        predictions, true_labels = [], []
        
        for batch in prediction_dataloader:
            batch = tuple(t.to(device) for t in batch)
            b_input_ids, b_input_mask, b_labels = batch
            
            with torch.no_grad():
                result = model(b_input_ids, 
                               token_type_ids=None, 
                               attention_mask=b_input_mask,
                               return_dict=True)
            
            logits = result.logits
            logits = logits.detach().cpu().numpy()
            label_ids = b_labels.to('cpu').numpy()
            
            predictions.append(logits)
            true_labels.append(label_ids)
        
        flat_predictions = np.concatenate(predictions, axis=0)
        flat_predictions = np.argmax(flat_predictions, axis=1).flatten()
        flat_true_labels = np.concatenate(true_labels, axis=0)
        
        f1 = f1_score(flat_true_labels, flat_predictions)
        print('F1 Score: %.3f' % f1)
        
        return f1
    
    else:
        # The original training and evaluation code
        data_df = pd.read_csv(tsv_path, sep='\t')
        data_df = data_df.fillna("")
        
        sen_w_feats = []
        labels = []
        
        for index, row in data_df.iterrows():
            fields = [str(row[field]) for field in combined_fields]
            combined = '[SEP]'.join(fields)
            sen_w_feats.append(combined)
            labels.append(row["Interact"])
        
        import numpy as np
        train_size = int(0.8 * len(data_df))
        val_size = int(0.1 * len(data_df))
        test_size = len(data_df) - (train_size + val_size)
        
        indices = np.arange(0, len(data_df))
        np.random.shuffle(indices)
        
        train_idx = indices[0:train_size]
        val_idx = indices[train_size:(train_size + val_size)]
        test_idx = indices[(train_size + val_size):]
        
        import random
        random.seed(42)
        np.random.seed(42)
        torch.manual_seed(42)
        torch.cuda.manual_seed_all(42)
        
        train_dataset = TensorDataset(input_ids[train_idx], attention_masks[train_idx], labels[train_idx])
        val_dataset = TensorDataset(input_ids[val_idx], attention_masks[val_idx], labels[val_idx])
        test_dataset = TensorDataset(input_ids[test_idx], attention_masks[test_idx], labels[test_idx])
        
        train_dataloader = DataLoader(train_dataset, sampler=RandomSampler(train_dataset), batch_size=8)
        validation_dataloader = DataLoader(val_dataset, sampler=SequentialSampler(val_dataset), batch_size=8)
        
        from transformers import AdamW, get_linear_schedule_with_warmup, BertTokenizer, BertForSequenceClassification
        
        tokenizer = BertTokenizer.from_pretrained('bert-base-uncased', do_lower_case=True)
        model = BertForSequenceClassification.from_pretrained("bert-base-uncased", num_labels=2)
        model.cuda()
        
        optimizer = AdamW(model.parameters(), lr=1e-5, eps=1e-8)
        total_steps = len(train_dataloader) * 10
        scheduler = get_linear_schedule_with_warmup(optimizer, num_warmup_steps=0, num_training_steps=total_steps)
        
        def flat_accuracy(preds, labels):
            pred_flat = np.argmax(preds, axis=1).flatten()
            labels_flat = labels.flatten()
            return np.sum(pred_flat == labels_flat) / len(labels_flat)
        
        import time
        import datetime
        import numpy as np
        
        def format_time(elapsed):
            elapsed_rounded = int(round((elapsed)))
            return str(datetime.timedelta(seconds=elapsed_rounded))
        
        training_stats = []
        total_t0 = time.time()
        
        for epoch_i in range(0, 10):
            print("")
            print('======== Epoch {:} / {:} ========'.format(epoch_i + 1, 10))
            print('Training...')
            
            t0 = time.time()
            total_train_loss = 0
            model.train()
            
            for step, batch in enumerate(train_dataloader):
                if step % 40 == 0 and not step == 0:
                    elapsed = format_time(time.time() - t0)
                    print('  Batch {:>5,}  of  {:>5,}.    Elapsed: {:}.'.format(step, len(train_dataloader), elapsed))
                
                b_input_ids = batch[0].to(device)
                b_input_mask = batch[1].to(device)
                b_labels = batch[2].to(device)
                
                model.zero_grad()
                result = model(b_input_ids, 
                               token_type_ids=None, 
                               attention_mask=b_input_mask, 
                               labels=b_labels,
                               return_dict=True)
                
                loss = result.loss
                logits = result.logits
                
                total_train_loss += loss.item()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()
                scheduler.step()
            
            avg_train_loss = total_train_loss / len(train_dataloader)
            training_time = format_time(time.time() - t0)
            
            print("")
            print("  Average training loss: {0:.2f}".format(avg_train_loss))
            print("  Training epoch took: {:}".format(training_time))
            
            print("")
            print("Running Validation...")
            
            t0 = time.time()
            model.eval()
            
            total_eval_accuracy = 0
            total_eval_loss = 0
            
            for batch in validation_dataloader:
                b_input_ids = batch[0].to(device)
                b_input_mask = batch[1].to(device)
                b_labels = batch[2].to(device)
                
                with torch.no_grad():
                    result = model(b_input_ids, 
                                   token_type_ids=None, 
                                   attention_mask=b_input_mask,
                                   labels=b_labels,
                                   return_dict=True)
                    
                loss = result.loss
                logits = result.logits
                
                total_eval_loss += loss.item()
                
                logits = logits.detach().cpu().numpy()
                label_ids = b_labels.to('cpu').numpy()
                
                total_eval_accuracy += flat_accuracy(logits, label_ids)
            
            avg_val_accuracy = total_eval_accuracy / len(validation_dataloader)
            print("  Accuracy: {0:.2f}".format(avg_val_accuracy))
            
            avg_val_loss = total_eval_loss / len(validation_dataloader)
            
            validation_time = format_time(time.time() - t0)
            
            print("  Validation Loss: {0:.2f}".format(avg_val_loss))
            print("  Validation took: {:}".format(validation_time))
            
            training_stats.append(
                {
                    'epoch': epoch_i + 1,
                    'Training Loss': avg_train_loss,
                    'Valid. Loss': avg_val_loss,
                    'Valid. Accur.': avg_val_accuracy,
                    'Training Time': training_time,
                    'Validation Time': validation_time
                }
            )
        
        print("")
        print("Training complete!")
        
        print("Total training took {:} (h:mm:ss)".format(format_time(time.time()-total_t0)))
        
        torch.save(model.state_dict(), model_name)
        
        df_stats = pd.DataFrame(data=training_stats)
        df_stats = df_stats.set_index('epoch')
        
        import matplotlib.pyplot as plt
        import seaborn as sns

        sns.set(style='darkgrid')
        sns.set(font_scale=1.5)
        plt.rcParams["figure.figsize"] = (12,6)

        plt.plot(df_stats['Training Loss'], 'b-o', label="Training")
        plt.plot(df_stats['Valid. Loss'], 'g-o', label="Validation")
        plt.title("Training & Validation Loss")
        plt.xlabel("Epoch")
        plt.ylabel("Loss")
        plt.legend()
        plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        plt.savefig(model_name.replace('.pth', '.png'))
        
        # Create a DataLoader to batch our test samples for us. We'll use a sequential
        # sampler this time--don't need this to be random!
        prediction_sampler = SequentialSampler(test_dataset)
        prediction_dataloader = DataLoader(test_dataset, sampler=prediction_sampler, batch_size=batch_size)

        print('Predicting labels for {:,} test sentences...'.format(len(test_dataset)))

        # Put model in evaluation mode
        model.eval()

        # Tracking variables 
        predictions , true_labels = [], []

        # Predict 
        for batch in prediction_dataloader:
            # Add batch to GPU
            batch = tuple(t.to(device) for t in batch)

            # Unpack the inputs from our dataloader
            b_input_ids, b_input_mask, b_labels = batch

            # Telling the model not to compute or store gradients, saving memory and 
            # speeding up prediction
            with torch.no_grad():        
                # Forward pass, calculate logit predictions.
                result = model(b_input_ids, 
                                token_type_ids=None, 
                                attention_mask=b_input_mask,
                                return_dict=True)

            logits = result.logits

            # Move logits and labels to CPU
            logits = logits.detach().cpu().numpy()
            label_ids = b_labels.to('cpu').numpy()

            # Store predictions and true labels
            predictions.append(logits)
            true_labels.append(label_ids)

        print('    DONE.')

        # Combine the results across all batches. 
        flat_predictions = np.concatenate(predictions, axis=0)

        # For each sample, pick the label (0 or 1) with the higher score.
        flat_predictions = np.argmax(flat_predictions, axis=1).flatten()

        # Combine the correct labels for each batch into a single list.
        flat_true_labels = np.concatenate(true_labels, axis=0)

        from sklearn.metrics import f1_score

        # Calculate the F1
        f1 = f1_score(flat_true_labels, flat_predictions)

        print('F1 Score: %.3f' % f1)
        df_stats['F1 Score'] = f1
        
        return df_stats
import sys
if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python train_bert.py <model_name> <tsv_path> <combined_field1> [<combined_field2> ...]")
        sys.exit(1)
    model_name = sys.argv[1]
    tsv_path = sys.argv[2]
    combined_fields = sys.argv[3:]
    main(model_name, tsv_path, combined_fields)
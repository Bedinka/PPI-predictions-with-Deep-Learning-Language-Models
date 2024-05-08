import json
import matplotlib.pyplot as plt
import numpy as np 

epochs_list = []
latent_dim_list = []
matrix_size_list = []
processed_sample_list = []
spearman_correlation_list = []
train_losses_list = []
val_losses_list = []
model_names = []

with open('summary_detailed.txt', 'r') as file:
    for line in file:
        if line.strip():  # Skip empty lines
            line = line.replace('nan', 'None')
            data = eval(line.strip())  
            model_names.append(data['model_name'])
            epochs_list.append(data['epochs'])
            latent_dim_list.append(data['latent_dim'])
            matrix_size_list.append(data['matrix_size'])
            processed_sample_list.append(data['processed_sample'])
            spearman_correlation_list.append(data['spearman_correlation'])
            train_losses_list.append(data['train_losses'])
            val_losses_list.append(data['val_losses'])

num_models = len(model_names)
colors = plt.cm.rainbow(np.linspace(0, 1, num_models))


# Plot spearman correlation vs epochs for all models
plt.figure(figsize=(16, 10))

for i, model_name in enumerate(model_names):
    for j in range(len(epochs_list)):
        plt.scatter(epochs_list[j], spearman_correlation_list[j], color=colors[j], marker='o', s=35, label=model_names[j])

plt.xlabel('Epochs')
plt.ylabel('Spearman Correlation')
plt.title('Spearman Correlation vs Epochs')
plt.grid(True)

# Create legend
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Model Name')
plt.subplots_adjust(right=0.6)

plt.tight_layout()
plt.savefig('spearman_correlation_vs_epochs.png')
plt.close()

# Plot spearman correlation vs matrix size for all models
plt.figure(figsize=(16, 10))
for i, model_name in enumerate(model_names):
    for j in range(len(matrix_size_list)):
        plt.scatter(matrix_size_list[j], spearman_correlation_list[j], color=colors[j], marker='o', s=35, label=model_names[j])

plt.xlabel('Matrix Size')
plt.ylabel('Spearman Correlation')
plt.title('Spearman Correlation vs Matrix Size')
plt.grid(True)

# Create legend
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Model Name')
plt.subplots_adjust(right=0.6)

plt.tight_layout()
plt.savefig('spearman_correlation_vs_matrix_size.png')
plt.close()

# Plot spearman correlation vs processed sample for all models
plt.figure(figsize=(16, 10))
for i, model_name in enumerate(model_names):
    for j in range(len(processed_sample_list)):
        plt.scatter(processed_sample_list[j], spearman_correlation_list[j], color=colors[j], marker='o', s=35, label=model_names[j])

plt.xlabel('Processed Sample')
plt.ylabel('Spearman Correlation')
plt.title('Spearman Correlation vs Processed Sample')
plt.grid(True)

# Create legend
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Model Name')
plt.subplots_adjust(right=0.6)

plt.tight_layout()
plt.savefig('spearman_correlation_vs_processed_sample.png')
plt.close()
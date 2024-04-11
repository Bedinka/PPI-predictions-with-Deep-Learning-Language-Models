import json
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np 


epochs_list = []
latent_dim_list = []
matrix_size_list = []
processed_sample_list = []
spearman_correlation_list = []
train_losses_list = []
val_losses_list = []
model_names = []
colors = cm.rainbow(np.linspace(0, 1, len(model_names)))

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

'''for index, model_data in enumerate(zip(epochs_list, latent_dim_list, matrix_size_list, processed_sample_list, train_losses_list, val_losses_list, model_names)):
    epochs, latent_dim, matrix_size, processed_sample, train_losses, val_losses, model_name = model_data
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, epochs + 1), train_losses, label='Train Loss', marker='o', markersize=3)
    plt.plot(range(1, epochs + 1), val_losses, label='Validation Loss', marker='o', markersize=3)
    plt.title(f'{model_name} Losses')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{model_name}_losses.png')
    plt.close()'''

for i, model_name in enumerate(model_names):
    # Plot each model's extracted data in a separate plot
    plt.figure(figsize=(10, 6))
    
    plt.subplot(1, 2, 1)
    plt.scatter(epochs_list[i], spearman_correlation_list[i], color=colors[i], marker='o', s=10, label=model_name)
    plt.xlabel('Epochs')
    plt.ylabel('Spearman Correlation')
    plt.title('Spearman Correlation vs Epochs')
    plt.grid(True)

    plt.subplot(1, 2, 2)
    for j in range(len(train_losses_list[i])):
        plt.plot(range(1, epochs_list[i] + 1), train_losses_list[i][j], color=colors[i], label=f'{model_name} Train Loss')
        plt.plot(range(1, epochs_list[i] + 1), val_losses_list[i][j], color=colors[i], label=f'{model_name} Validation Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.title(f'{model_name} Losses')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(f'combined_plot_{model_name}.png')
    plt.close()
    
# Plot spearman correlation, train losses, and val losses against each parameter for each model
plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
for i, model_name in enumerate(model_names):
    plt.scatter(epochs_list[i], spearman_correlation_list[i], color='blue', marker='o', s=10, label=model_name)
plt.xlabel('Epochs')
plt.ylabel('Spearman Correlation')
plt.title('Spearman Correlation vs Epochs')
plt.grid(True)

plt.subplot(2, 2, 2)
for i, model_name in enumerate(model_names):
    for j in range(len(train_losses_list[i])):
        plt.scatter(latent_dim_list[i], train_losses_list[i][j], color='red', label=f'{model_name} Train Loss', marker='o', s=10)
        plt.scatter(latent_dim_list[i], val_losses_list[i][j], color='green', label=f'{model_name} Validation Loss', marker='o', s=10)
plt.xlabel('Latent Dimension')
plt.ylabel('Loss')
plt.title('Losses vs Latent Dimension')
plt.grid(True)

plt.subplot(2, 2, 3)
for i, model_name in enumerate(model_names):
    plt.scatter(matrix_size_list[i], spearman_correlation_list[i], color='blue', marker='o', s=10, label=model_name)
plt.xlabel('Matrix Size')
plt.ylabel('Spearman Correlation')
plt.title('Spearman Correlation vs Matrix Size')
plt.grid(True)

plt.subplot(2, 2, 4)
for i, model_name in enumerate(model_names):
    plt.scatter(processed_sample_list[i], spearman_correlation_list[i], color='blue', marker='o', s=10, label=model_name)
plt.xlabel('Processed Sample')
plt.ylabel('Spearman Correlation')
plt.title('Spearman Correlation vs Processed Sample')
plt.grid(True)

plt.tight_layout()
plt.savefig('combined_plot.png')
plt.close()
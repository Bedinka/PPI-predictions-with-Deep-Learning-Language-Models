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

"""for index, model_data in enumerate(zip(epochs_list, latent_dim_list, matrix_size_list, processed_sample_list, train_losses_list, val_losses_list, model_names)):
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
    plt.close()"""


num_models = len(model_names)
colors = plt.cm.rainbow(np.linspace(0, 1, num_models))

'''for i, model_name in enumerate(model_names):
    fig, ax = plt.subplots(figsize=(16,12))
    fig.subplots_adjust(right=0.6)

    for j in range(len(epochs_list)):
        ax.scatter(epochs_list[j], spearman_correlation_list[j], color=colors[j], marker='o', s=35, label=model_names[j])
        
    ax.set_xlabel('Epochs')
    ax.set_ylabel('Spearman Correlation')
    ax.set_title('Spearman Correlation vs Epochs')
    ax.grid(True)
    
    # Create legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1.05, 0.5), title='Model Name')

    # Save the plot
    plt.savefig(f'{model_name}_plot.png')
    plt.close()'''
# Plot spearman correlation, train losses, and val losses against each parameter for each model
plt.figure(figsize=(17, 13))
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.4, hspace=0.4)
plt.subplot(4, 1, 1)
for i, model_name in enumerate(model_names):
    plt.scatter(epochs_list[i], spearman_correlation_list[i], color=colors[i], marker='o', s=30, label=model_name)
plt.xlabel('Epochs')
plt.ylabel('Spearman Correlation')
plt.title('Spearman Correlation vs Epochs')
plt.grid(True)

plt.subplot(4, 1, 2)
for i, model_name in enumerate(model_names):
    for j in range(len(train_losses_list[i])):
        plt.scatter(latent_dim_list[i], train_losses_list[i][j], color='red', label=f'{model_name} Train Loss', marker='o', s=20)
        plt.scatter(latent_dim_list[i], val_losses_list[i][j], color='green', label=f'{model_name} Validation Loss', marker='o', s=20)
plt.xlabel('Latent Dimension')
plt.ylabel('Loss')
plt.title('Losses vs Latent Dimension')
plt.grid(True)

plt.subplot(4, 1, 3)
for i, model_name in enumerate(model_names):
    plt.scatter(matrix_size_list[i], spearman_correlation_list[i], color=colors[i], marker='o', s=30, label=model_name)
plt.xlabel('Matrix Size')
plt.ylabel('Spearman Correlation')
plt.title('Spearman Correlation vs Matrix Size')
plt.grid(True)

plt.subplot(4, 1, 4)
for i, model_name in enumerate(model_names):
    plt.scatter(processed_sample_list[i], spearman_correlation_list[i], color=colors[i], marker='o', s=30, label=model_name)
plt.xlabel('Processed Sample')
plt.ylabel('Spearman Correlation')
plt.title('Spearman Correlation vs Processed Sample')
plt.grid(True)

plt.legend(loc= 'lower left', bbox_to_anchor=(1 ,0.5), title='Model Name', fontsize= 'smaller')
plt.subplots_adjust(right=0.7)
plt.tight_layout()
plt.savefig('combined_plot.png')
plt.close()
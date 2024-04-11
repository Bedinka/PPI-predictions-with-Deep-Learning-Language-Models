import json
import matplotlib.pyplot as plt


with open('summary.txt', 'r') as file:
    lines = file.readlines()
val_losses_values = []
spearman_correlation_values = []
spearman_correlation_dim1_2_values = []
spearman_p_value_dim1_2_values = []
model_names = []

for line in lines:

    model_name, data_str = line.split("\n")[0].split(" ", 1)
    
    data = json.loads(data_str)
    
    model_names.append(model_name)
    val_losses_values.append(data['val_losses'])
    spearman_correlation_values.append(data['spearman_correlation'])
    spearman_correlation_dim1_2_values.append(data['spearman_correlation_dim1_2'])
    spearman_p_value_dim1_2_values.append(data['spearman_p_value_dim1_2'])

plt.figure(figsize=(10, 6))

for i, val_losses in enumerate(val_losses_values):
    plt.plot(val_losses, label=model_names[i])

plt.title('Val Losses')
plt.legend()
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.grid(True)

plt.tight_layout()

plt.savefig('val_losses_plot.png')

plt.close()

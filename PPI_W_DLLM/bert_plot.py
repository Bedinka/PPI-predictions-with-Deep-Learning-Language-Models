import matplotlib.pyplot as plt

# Data from the BERT model testing
model_names = [
    "bert_attrnum_1",
    "bert_attrnum_2",
    "bert_attrnum_3",
    "bert_attrnum_4"
]

accuracies = [0.6150, 0.6429, 0.6265, 0.6435]
sensitivities = [0.6890, 0.8110, 0.8171, 0.7866]
specificities = [0.6069, 0.6243, 0.6055, 0.6277]
auc_values = [0.6487, 0.7831, 0.7730, 0.7759]

# Plotting
plt.figure(figsize=(12, 8))

# Plot Accuracy
plt.subplot(2, 2, 1)
plt.plot(model_names, accuracies, marker='o', linestyle='-', color='b')
plt.title('Accuracy')
plt.xlabel('Model')
plt.ylabel('Accuracy')
plt.xticks(rotation=45)
plt.grid(True)

# Plot Sensitivity
plt.subplot(2, 2, 2)
plt.plot(model_names, sensitivities, marker='o', linestyle='-', color='g')
plt.title('Sensitivity')
plt.xlabel('Model')
plt.ylabel('Sensitivity')
plt.xticks(rotation=45)
plt.grid(True)

# Plot Specificity
plt.subplot(2, 2, 3)
plt.plot(model_names, specificities, marker='o', linestyle='-', color='r')
plt.title('Specificity')
plt.xlabel('Model')
plt.ylabel('Specificity')
plt.xticks(rotation=45)
plt.grid(True)

# Plot AUC
plt.subplot(2, 2, 4)
plt.plot(model_names, auc_values, marker='o', linestyle='-', color='m')
plt.title('AUC')
plt.xlabel('Model')
plt.ylabel('AUC')
plt.xticks(rotation=45)
plt.grid(True)

plt.tight_layout()
plt.show()

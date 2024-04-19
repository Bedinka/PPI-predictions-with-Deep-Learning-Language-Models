import seaborn as sns
import matplotlib.pyplot as plt

import data_collect

df = data_collect.main()
corr = df.corr()
plt.figure(figsize=(10, 8))
sns.heatmap(corr, annot=True, cmap='coolwarm', fmt=".2f", linewidths=.5)
plt.title('Correlation Heatmap')
plt.show()
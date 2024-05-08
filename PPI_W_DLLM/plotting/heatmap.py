import YangLabIntern.PPI_W_DLLM.Autoencoder.autoencoder_data_collect as autoencoder_data_collect
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#df = data_collect.main()
df = pd.read_csv("results.csv")
df_numeric = df.drop(columns=["Model"])
corr = df_numeric.corr()

plt.figure(figsize=(10, 8))
sns.heatmap(corr, annot=True, cmap='coolwarm', fmt=".2f", linewidths=.5)
plt.title('Correlation Heatmap')
plt.show()

plt.figure(figsize=(12, 8))
sns.scatterplot(data=df, x="Latent Dimension", y="Matrix Size", hue="Spearman Correlation", size="Processed Sample", sizes=(50, 200), palette="coolwarm", legend="full")
plt.title("Correlation of Keras Models")
plt.xlabel("Latent Dimension")
plt.ylabel("Matrix Size")
plt.grid(True)
plt.legend(title="Spearman Correlation")
plt.show()
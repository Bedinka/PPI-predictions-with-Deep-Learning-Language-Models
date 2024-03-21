import feature_extraction
import numpy as np
import pandas as pd
import tensorflow as tf

from sklearn.metrics import accuracy_score, precision_score, recall_score
from sklearn.model_selection import train_test_split
from tensorflow.keras import layers, losses
from tensorflow.keras.datasets import fashion_mnist
from tensorflow.keras.models import Model

interacting_prot = []
interacting_prot = feature_extraction.main()

for protein in interacting_prot:
    print(protein.distance_matrices_CA)
    dist_matrix = protein.distance_matrices_CA
    num_rows = len(dist_matrix)
    mid_point = num_rows // 2
    dist_ca_train = dist_matrix[:mid_point]
    dist_ca_test = dist_matrix[mid_point:]
    """
    (dist_ca_ab_train, _)= (dist_ca_ab_test, _)= protein.distance_matrices_CA_AB
    (dist_ca_ab_sub_train, _)= (dist_ca_ab_sub_test, _)= protein.ca_submatrices
    (dist_mean_train, _)= (dist_mean_test, _)= protein.distance_matrices_mean
    (dist_mean_ab_train, _)= (dist_mean_ab_test, _)= protein.distance_matrices_mean_AB
    (dist_mean_ab_sub_train, _)= (dist_mean_ab_sub_test, _)= protein.mean_submatrices
    """
dist_ca_train = dist_ca_train.astype('float32') / 255.
dist_ca_test = dist_ca_test.astype('float32') / 255.

print (dist_ca_train.shape)
print (dist_ca_test.shape)

class Autoencoder(Model):

  def __init__(self, latent_dim, shape):
    super(Autoencoder, self).__init__()
    self.latent_dim = latent_dim
    self.shape = shape
    self.encoder = tf.keras.Sequential([
      layers.Flatten(),
      layers.Dense(latent_dim, activation='relu'),
    ])
    self.decoder = tf.keras.Sequential([
      layers.Dense(tf.math.reduce_prod(shape), activation='sigmoid'),
      layers.Reshape(shape)
    ])

  def call(self, x):
    encoded = self.encoder(x)
    decoded = self.decoder(encoded)
    return decoded


shape = dist_ca_train.shape[1:]
latent_dim = 100
autoencoder = Autoencoder(latent_dim, shape)

autoencoder.compile(optimizer='adam', loss=losses.MeanSquaredError())

autoencoder.fit(dist_ca_train, dist_ca_train,
                epochs=10,
                shuffle=True,
                validation_data=(dist_ca_test, dist_ca_test))


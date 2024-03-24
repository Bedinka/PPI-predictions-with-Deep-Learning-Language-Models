import feature_extraction
import numpy as np
import pandas as pd
import tensorflow as tf
import os

from sklearn.metrics import accuracy_score, precision_score, recall_score
from sklearn.model_selection import train_test_split
from tensorflow.keras import layers, losses
from tensorflow.keras.datasets import fashion_mnist
from tensorflow.keras.models import Model

os.environ['CUDA_VISIBLE_DEVICES'] = "0"


interacting_prot = []
interacting_prot = feature_extraction.main()


def mean_matrix_vec(proteins):
  for protein in proteins:
    if protein.chainID == 'A':
      #sub_matrix = np.array([protein.mean_submatrices, protein.ca_submatrices])
      sub_values_a = np.array([protein.mean_submatrices],dtype=np.float32)
      print(sub_values_a.shape)
      re_a = sub_values_a.reshape(-1, 7, 7)
    else: 
      sub_values_b = np.array([protein.mean_submatrices], dtype=np.float32)
      print(sub_values_b.shape)
      re_b= sub_values_b.reshape(-1, 7, 7)
  return re_a, re_b 

def ca_matrix_vec(proteins):
  for protein in proteins:
    if protein.chainID == 'A':
      sub_values_a = np.array([protein.ca_submatrices], dtype=np.float32)
      print(sub_values_a.shape)
      re_a = sub_values_a.reshape(-1, 7, 7)
    else: 
      sub_values_b = np.array([protein.ca_submatrices], dtype=np.float32)
      print(sub_values_b.shape)
      re_b= sub_values_b.reshape(-1, 7, 7)
  return re_a,re_b 


re_a, re_b = mean_matrix_vec(interacting_prot)
ca_a, ca_b = ca_matrix_vec(interacting_prot)
final_train_vec =  np.concatenate((re_a, ca_a), axis=0)
final_test_vec = np.concatenate((re_b, ca_b), axis=0) 
    #sub_matrix_2 = np. array ([[protein.mean_submatrices],[protein.ca_submatrices]], dtype = float )

dist_ca_train =  final_train_vec
dist_ca_test = final_test_vec
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
latent_dim = 10
autoencoder = Autoencoder(latent_dim, shape)

autoencoder.compile(optimizer='adam', loss=losses.MeanSquaredError())

autoencoder.fit(dist_ca_train, dist_ca_train,
                epochs=10,
                shuffle=True,
                validation_data=(dist_ca_test, dist_ca_test))


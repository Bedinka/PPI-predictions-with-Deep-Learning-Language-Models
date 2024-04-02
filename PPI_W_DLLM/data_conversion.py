import feature_extraction
import numpy as np
import pandas as pd
import tensorflow as tf
import os
import torch 

from sklearn.metrics import accuracy_score, precision_score, recall_score
from sklearn.model_selection import train_test_split
from tensorflow.keras import layers, losses
from tensorflow.keras.models import Model

from scipy import stats
from random import randint
#os.environ['CUDA_VISIBLE_DEVICES'] = "0"

interacting_prot = []
interacting_prot = feature_extraction.main()

def mean_matrix_vec(interacting_prot): 
  re_a = []
  re_b = []
  for protein in interacting_prot:
    if protein.chainID == 'A':
      sub_values_a = np.array([protein.mean_submatrices])
      print(sub_values_a.shape)
      tensor1 = torch.tensor(sub_values_a.reshape(-1, 7, 7))
      re_a.append(sub_values_a.reshape(-1, 7, 7))
    else: 
      sub_values_b = np.array([protein.mean_submatrices])
      print(sub_values_b.shape)
      tensor2 = torch.tensor(sub_values_b.reshape(-1, 7, 7))
      re_b.append(sub_values_b.reshape(-1, 7, 7))
  
  return re_a, re_b , tensor1 , tensor2
    
def ca_matrix_vec(proteins):
  for protein in proteins:
    if protein.chainID == 'A':
      sub_values_a = np.array([protein.ca_submatrices])
      print(sub_values_a.shape)
      ca_a = sub_values_a.reshape(-1, 7, 7)
    else: 
      sub_values_b = np.array([protein.ca_submatrices])
      print(sub_values_b.shape)
      ca_b= sub_values_b.reshape(-1, 7, 7)
  return ca_a,ca_b 


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
      layers.Dense(tf.math.reduce_prod(shape).numpy(), activation='sigmoid'),
      layers.Reshape(shape)
    ])

  def call(self, x):
    encoded = self.encoder(x)
    decoded = self.decoder(encoded)
    return decoded


def main():

  re_a, re_b , tensor1 , tensor2 = mean_matrix_vec(interacting_prot)
  ca_a, ca_b =ca_matrix_vec(interacting_prot)
  finals_train_vec = np.concatenate(re_a, axis=0)
  finals_test_vec = np.concatenate(re_b, axis=0) 
  dist_ca_train =  finals_train_vec
  dist_ca_test =  finals_test_vec
  dist_ca_train = dist_ca_train.astype('float32') / 255.
  dist_ca_test = dist_ca_test.astype('float32') / 255.
  print (dist_ca_train.shape)
  print (dist_ca_test.shape)
  shape = dist_ca_train.shape[1:]
  latent_dim = 1
  autoencoder = Autoencoder(latent_dim, shape)

  autoencoder.compile(optimizer='adam', loss=losses.MeanSquaredError())

  autoencoder.fit(dist_ca_train, dist_ca_train,
                  epochs=10,
                  shuffle=True,
                  validation_data=(dist_ca_test, dist_ca_test))
  encoded_vectors_train = autoencoder.encoder(dist_ca_train)
  encoded_vectors_train_2d = dist_ca_train.reshape(dist_ca_train.shape[0], -1)
  tensor1_2d = tensor1.reshape(tensor1.shape[0], -1)
  correlation_coefficient = np.corrcoef(encoded_vectors_train_2d, tensor1_2d)[0, 1]
  '''
  1 indicates a perfect positive linear relationship.
  -1 indicates a perfect negative linear relationship.
  0 indicates no linear relationship.'''
  print(f"Correlation: {correlation_coefficient}")
  print(encoded_vectors_train.shape)
  print(encoded_vectors_train)
  encoded_vectors_test = autoencoder.encoder(dist_ca_test)
  encoded_vectors_train_tensor = tf.constant(encoded_vectors_train)
  encoded_vectors_test_tensor = tf.constant(encoded_vectors_test)

  print(dist_ca_train.shape)
  print(encoded_vectors_train.shape)
  print(dist_ca_train[1,].shape)
  print(dist_ca_train[2,].shape)

  X = []
  Y = []
  
  for i in range(100):
    j = 100 #randint(1, 1000)

    dist1 = np.sum(np.abs(dist_ca_train[i,] - dist_ca_train[i+j,]))
    #dist2 = np.sum(np.abs(encoded_vectors_train[i,] - encoded_vectors_train[i+j,]))

    #dist1 = euclidean_distance(dist_ca_train[i,], dist_ca_train[i+1,])
    X.append(dist1)
    #Y.append(dist2)

    #print('dist1', dist1)
    #print('dist2', dist2)
    
    
  print(encoded_vectors_train[1,].shape)
  print(encoded_vectors_train[2,].shape)
  
  for i in range(100):
    j = 100#randint(1, 1000)
    
    #dist1 = np.sum(np.abs(dist_ca_train[i,] - dist_ca_train[i+j,]))
    dist2 = np.sum(np.abs(encoded_vectors_train[i,] - encoded_vectors_train[i+j,]))

    #dist1 = euclidean_distance(dist_ca_train[i,], dist_ca_train[i+1,])
    #X.append(dist1)
    Y.append(dist2)

    #print('dist1', dist1)
    #print('dist2', dist2)

  
  print(stats.spearmanr(X, Y))
  autoencoder.save(f'autoencoder_models/dina_model_dim{latent_dim}.keras')
  return encoded_vectors_train_tensor


if __name__ == "__main__":
    main()
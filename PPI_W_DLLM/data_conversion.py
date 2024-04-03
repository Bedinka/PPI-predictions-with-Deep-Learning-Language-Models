import feature_extraction
import numpy as np
import pandas as pd
import tensorflow as tf
import os
import torch 
from scipy import stats
import matplotlib.pyplot as plt

from sklearn.metrics import accuracy_score, precision_score, recall_score
from sklearn.model_selection import train_test_split
from tensorflow.keras import layers, losses
from tensorflow.keras.models import Model
from keras.callbacks import Callback

#from autoencoder import Autoencoder

from scipy import stats
from random import randint
#os.environ['CUDA_VISIBLE_DEVICES'] = "0"

interacting_prot = []
interacting_prot = feature_extraction.main()

@tf.keras.utils.register_keras_serializable()
class Autoencoder(Model):

    def __init__(self, latent_dim, shape, **kwargs):
        super(Autoencoder, self).__init__(**kwargs)
        self.latent_dim = latent_dim
        self.shape = shape
        self.encoder = tf.keras.Sequential([
            tf.keras.layers.Flatten(),
            tf.keras.layers.Dense(latent_dim, activation='relu'),
        ])
        self.decoder = tf.keras.Sequential([
            tf.keras.layers.Dense(tf.math.reduce_prod(shape).numpy(), activation='sigmoid'),
            tf.keras.layers.Reshape(shape)
        ])

    def call(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded
    def get_config(self):
      config = super().get_config()
      config.update(
          {
            "latent_dim": self.latent_dim,
            "shape": self.shape
          }
      )
      return config
    
    def custom_objects():
    
      return {"Autoencoder": Autoencoder}
    
def mean_matrix_vec(interacting_prot): 
  re = []
  re_b = []
  sub_index = []
  for protein in interacting_prot:
    if protein.chainID == 'A':
      sub_values_a = np.array([protein.mean_submatrices])
      print(sub_values_a.shape)
      tensor1 = torch.tensor(sub_values_a.reshape(-1, 7, 7))
      re.append(sub_values_a.reshape(-1, 7, 7))
      sub_index.append(np.array([protein.sub_res_index]))
    else: 
      sub_values_b = np.array([protein.mean_submatrices])
      print(sub_values_b.shape)
      tensor2 = torch.tensor(sub_values_b.reshape(-1, 7, 7))
      re.append(sub_values_b.reshape(-1, 7, 7))
      sub_index.append(np.array([protein.sub_res_index]))
  return re
    
def ca_matrix_vec(proteins):
  ca = []
  ca_b =[]
  for protein in proteins:
    if protein.chainID == 'A':
      sub_values_a = np.array([protein.ca_submatrices])
      print(sub_values_a.shape)
      ca.append(sub_values_a.reshape(-1, 7, 7))
    else: 
      sub_values_b = np.array([protein.ca_submatrices])
      print(sub_values_b.shape)
      ca.append(sub_values_b.reshape(-1, 7, 7))
  return ca


def spearman(dist_ca_train,encoded_vectors_train ):
  X = []
  Y = []
  random  = randint(1, 1000)
  for i in range(1000):
    j = 100
    dist1 = np.sum(np.abs(dist_ca_train[i,] - dist_ca_train[i+j,]))
    X.append(dist1)

  print(encoded_vectors_train[1,].shape)
  print(encoded_vectors_train[2,].shape)
  
  for i in range(1000):
    j = 100
    dist2 = np.sum(np.abs(encoded_vectors_train[i,] - encoded_vectors_train[i+j,]))
    Y.append(dist2)
  
  print(stats.spearmanr(X, Y))
  return X, Y 

def plot(m1, m2):
  m1 = np.array(m1)
  m2 = np.array(m2)
  print(m1.shape[0])
  print(m2.shape[0])
  xmin = m1.min()
  xmax = m1.max()
  ymin = m2.min()
  ymax = m2.max()
  X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
  positions = np.vstack([X.ravel(), Y.ravel()])
  values = np.vstack([m1, m2])
  kernel = stats.gaussian_kde(values)
  Z = np.reshape(kernel(positions).T, X.shape)
                
  fig, ax = plt.subplots()
  ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,extent=[xmin, xmax, ymin, ymax])
  ax.plot(m1, m2, 'k.', markersize=2)
  ax.set_xlim([xmin, xmax])
  ax.set_ylim([ymin, ymax])
 
  ax.set_xlim([xmin, xmax])
  ax.set_ylim([ymin, ymax])
  plt.show()


def main():
  train_losses = []
  val_losses = []

  re_a = mean_matrix_vec(interacting_prot)
  ca_a =ca_matrix_vec(interacting_prot)
  finals_train_vec = np.concatenate(re_a, axis=0)
  finals_test_vec = np.concatenate(ca_a, axis=0) 
  dist_ca_train =  finals_train_vec
  dist_ca_test =  finals_test_vec
  dist_ca_train = dist_ca_train.astype('float32') / 255.
  dist_ca_test = dist_ca_test.astype('float32') / 255.
  print (dist_ca_train.shape)
  print (dist_ca_test.shape)
  shape = dist_ca_train.shape[1:]
  latent_dim = 2
  autoencoder = tf.keras.models.load_model('dina_model.keras', custom_objects={"Autoencoder": Autoencoder} ) #Autoencoder(latent_dim, shape)
  #tf.keras.models.load_model('/home/pc550/Documents/PPI_W_DLLM/autoencoder_models/dina_model_dim2.keras') 
  autoencoder.compile(optimizer='adam', loss=losses.MeanSquaredError())
  logs = Callback()

  autoencoder.fit(dist_ca_train, dist_ca_train,
                  epochs=10,
                  shuffle=True,
                  validation_data=(dist_ca_test, dist_ca_test),
                  callbacks=[logs])
  encoded_vectors_train = autoencoder.encoder(dist_ca_train)
  print(encoded_vectors_train.shape)
  print(encoded_vectors_train)
  encoded_vectors_test = autoencoder.encoder(dist_ca_test)
  encoded_vectors_train_tensor = tf.constant(encoded_vectors_train)
  encoded_vectors_test_tensor = tf.constant(encoded_vectors_test)

  print(dist_ca_train.shape)
  print(encoded_vectors_train.shape)
  print(dist_ca_train[1,].shape)
  print(dist_ca_train[2,].shape)

  x, y = spearman ( dist_ca_train, encoded_vectors_train)
  plot(x , y )
  #autoencoder.save('dina_model.keras')
  return encoded_vectors_train_tensor

if __name__ == "__main__":
    main()
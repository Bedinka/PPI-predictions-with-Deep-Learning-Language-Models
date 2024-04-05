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
from tensorflow.keras import layers, losses, callbacks
from tensorflow.keras.models import Model
from keras.callbacks import Callback
from sklearn.decomposition import PCA


#from autoencoder import Autoencoder

from scipy import stats
from random import randint
#os.environ['CUDA_VISIBLE_DEVICES'] = "0"


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
    
class LossHistory(callbacks.Callback):
    def __init__(self):
        super().__init__()
        self.train_losses = []
        self.val_losses = []

    def on_epoch_end(self, epoch, logs=None):
        self.train_losses.append(logs.get('loss'))
        self.val_losses.append(logs.get('val_loss'))
   
def mean_matrix_vec(interacting_prot, size): 
  re = []
  sub_index = []
  for protein in interacting_prot:
    sub_values= np.array([protein.mean_submatrices])
    print(sub_values.shape)
    re.extend(sub_values.reshape(-1, size+1, size+1))
    #sub_index.append(np.array([protein.sub_res_index]))
    
  return re
def ca_matrix_vec(proteins, size):
  ca = []
  for protein in proteins:
      sub_values = np.array([protein.ca_submatrices])
      print(sub_values.shape)
      ca.extend(sub_values.reshape(-1, size+1, size+1))
  return ca


def spearman(dist_ca_train,encoded_vectors_train, ranges , int  ):
  X = []
  Y = []
  random  = randint(1, 1000)
  for i in range(ranges):
    j = int
    dist1 = np.sum(np.abs(dist_ca_train[i,] - dist_ca_train[i+j,]))
    X.append(dist1)
  
  for i in range(ranges):
    j = int
    dist2 = np.sum(np.abs(encoded_vectors_train[i,] - encoded_vectors_train[i+j,]))
    Y.append(dist2)
  
  correlation, p_value =stats.spearmanr(X, Y)
  print(encoded_vectors_train[1,])

  for i in range(5):  
    print(encoded_vectors_train[i, 0])  
    print(encoded_vectors_train[i, 1])

  return X, Y , correlation, p_value 

def plot(encoded):
  
  m1 = np.array(encoded[1])
  m2 = np.array(encoded[0])
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
  

def main(latent_dim, model_path, processed_sample, size, epochs=10):
  
  import feature_extraction
  interacting_prot = feature_extraction.main(processed_sample, size)
  ranges= 2000
  int = 100
  loss_history = LossHistory()
  re_a = mean_matrix_vec(interacting_prot, size)
  ca_a =ca_matrix_vec(interacting_prot, size)
  dist_ca_train = np.array(re_a)
  dist_ca_test =  np.array(ca_a)
  #dist_ca_train = np.concatenate(re_a, axis=0)
  #dist_ca_test =  np.concatenate(ca_a, axis=0)
  dist_ca_train = dist_ca_train.astype('float32') / 255.
  dist_ca_test = dist_ca_test.astype('float32') / 255.
  print (dist_ca_train.shape)
  print (dist_ca_test.shape)
  shape = dist_ca_train.shape[1:]
  autoencoder = Autoencoder(latent_dim, shape)
  #tf.keras.models.load_model(model_path, custom_objects={"Autoencoder": Autoencoder} ) 
  autoencoder.compile(optimizer='adam', loss=losses.MeanSquaredError())
  autoencoder.fit(dist_ca_train, dist_ca_train,
                  epochs=10,
                  shuffle=True,
                  validation_data=(dist_ca_test, dist_ca_test),
                  callbacks=[loss_history])
  encoded_vectors_train = autoencoder.encoder(dist_ca_train)
  print(encoded_vectors_train.shape)
  print(encoded_vectors_train)
  encoded_vectors_test = autoencoder.encoder(dist_ca_test)

  print(dist_ca_train.shape)
  print(encoded_vectors_train.shape)

  x, y , correlation, p_value = spearman ( dist_ca_train, encoded_vectors_train , ranges , int )
  
  plot(encoded_vectors_train )
  autoencoder.save('dina_model_1.keras')
  
  collected_data = {
     "train_losses": loss_history.train_losses,
      "val_losses": loss_history.val_losses,
      "spearman_correlation": correlation ,
      "spearman_p_value": p_value
  }
  

  
  return collected_data

if __name__ == "__main__":
    latent_dim = 10
    model_path = 'dina_model.keras'
    processed_sample = 30
    size = 7

    main(latent_dim ,model_path, processed_sample ,size ) #latent_dim ,model_path, processed_sample ,size
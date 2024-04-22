import feature_extraction
import numpy as np
import pandas as pd
import tensorflow as tf
import os
import torch 
from scipy import stats
import matplotlib.pyplot as plt
from sklearn.preprocessing import OneHotEncoder
import blosum as bl

from sklearn.metrics import accuracy_score, precision_score, recall_score
from sklearn.model_selection import train_test_split
from tensorflow.keras import layers, losses, callbacks
from tensorflow.keras.models import Model
from keras.callbacks import Callback
from sklearn.decomposition import PCA

from scipy import stats
from random import randint


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
  for protein in interacting_prot:
    sub_values= np.array([protein.mean_submatrices])
    #print(sub_values.shape)
    re.extend(sub_values.reshape(-1, size, size))
    #sub_index.append(np.array([protein.sub_res_index]))
    
  return re

def ca_matrix_vec(proteins, size):
  ca = []
  for protein in proteins:
      sub_values = np.array([protein.ca_submatrices])
      #print(sub_values.shape)
      ca.extend(sub_values.reshape(-1, size, size))
  return ca

def one_hot_encoding(submatrix):
  encoder = OneHotEncoder()
  one_hot_encoded = encoder.fit_transform(submatrix).toarray()
  print(submatrix)
  print(one_hot_encoded)
  np.savetxt('onehot.txt', one_hot_encoded)
  return one_hot_encoded

def blosum():
  matrix = bl.BLOSUM(62)
  val = matrix["A"]["Y"]

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
  #print( "Correlation between distances of encoded vectors and distance matrices", correlation, p_value )  
  
  #print(encoded_vectors_train[1,])
  #print(encoded_vectors_train[1,].shape)
  #print(stats.spearmanr(encoded_vectors_train[0,:],encoded_vectors_train[1,:]))

  try:
    DIM1 = []
    DIM2 = []
    for i in range(10000):  
      #print(encoded_vectors_train[i, 0])  
      #print(encoded_vectors_train[i, 1])
      DIM1.append(encoded_vectors_train[i, 0])
      DIM2.append(encoded_vectors_train[i, 1])
    
    correlation2, p_value2 =stats.spearmanr(DIM1, DIM2)
  except:
    correlation2 = None
    p_value2 = None

  #print( correlation2, p_value2 )  
  return X, Y , correlation, p_value, correlation2, p_value2 

def plot(encoded, model_name):
  try:
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
    plt.savefig(f'{model_name}_plot.png')
  except np.linalg.LinAlgError as e:
        print("Error: ", e)
  

def main(latent_dim, model_name, processed_sample, size, SAVE, epoch):
  
  import feature_extraction
  interacting_prot = feature_extraction.main(processed_sample, size)
  ranges= 2000
  int = 500
  loss_history = LossHistory()
  re_a = mean_matrix_vec(interacting_prot, size)
  ca_a =ca_matrix_vec(interacting_prot, size)
  
  dist_ca_train = np.array(re_a)
  dist_ca_test =  np.array(ca_a)
  dist_ca_train = dist_ca_train.astype('float32') / 255.
  dist_ca_test = dist_ca_test.astype('float32') / 255.
  #print (dist_ca_train.shape)
  #print (dist_ca_test.shape)
  shape = dist_ca_train.shape[1:]
  if os.path.exists(model_name):
    autoencoder = tf.keras.models.load_model(model_name, custom_objects={"Autoencoder": Autoencoder} )
  else: 
    autoencoder = Autoencoder(latent_dim, shape)
    autoencoder.compile(optimizer='adam', loss=losses.MeanSquaredError())
    autoencoder.fit(dist_ca_train, dist_ca_train,
                    epochs=epoch,
                    shuffle=True,
                    validation_data=(dist_ca_test, dist_ca_test),
                    callbacks=[loss_history])

  """if SAVE:
    autoencoder = Autoencoder(latent_dim, shape)
    autoencoder.compile(optimizer='adam', loss=losses.MeanSquaredError())
    autoencoder.fit(dist_ca_train, dist_ca_train,
                    epochs=epoch,
                    shuffle=True,
                    validation_data=(dist_ca_test, dist_ca_test),
                    callbacks=[loss_history])
  else:  
    autoencoder = tf.keras.models.load_model(model_name, custom_objects={"Autoencoder": Autoencoder} ) """
  
  encoded_vectors_train = autoencoder.encoder(dist_ca_train).numpy()
  #print(encoded_vectors_train.shape)
  #print(encoded_vectors_train)
  encoded_vectors_test = autoencoder.encoder(dist_ca_test).numpy()

  #print(dist_ca_train.shape)
  #print(encoded_vectors_train.shape)

  x, y , correlation, p_value, correlation2, p_value2  = spearman ( dist_ca_train, encoded_vectors_train , ranges , int )
  x_test, y_test , correlation_test, p_value_test, correlation2_test, p_value2_test  = spearman ( dist_ca_test, encoded_vectors_test , ranges , int )
  
  print("#################", model_name)
  print(correlation, p_value, correlation2, p_value2 )
  print(correlation_test, p_value_test, correlation2_test, p_value2_test )
  #one_hot_train = one_hot_encoding(encoded_vectors_train)
  plot(encoded_vectors_train, model_name)
  
  if SAVE:
    tf.saved_model.save(autoencoder, model_name)
  
  collected_data = {
     "model_name": model_name,
     "epochs" : epoch,
     "latent_dim" : latent_dim,
     "matrix_size" : size,
     "processed_sample" : processed_sample,
     "train_losses": loss_history.train_losses,
      "val_losses": loss_history.val_losses,
      "spearman_correlation": correlation ,
      "spearman_p_value": p_value,
      "spearman_correlation_dim1_2": correlation2,
      "spearman_p_value_dim1_2": p_value2
    }
  
  return encoded_vectors_train #, collected_data

if __name__ == "__main__":
    '''latent_dim = 
    model_name = 'dina_model_js.keras'
    processed_sample = 1
    size = 10'''
    #SAVE = True

    '''fout = open("summary.txt", "a")
    for i in range(10):
      model_path = 'dina_model_sample_%d_dim_%d_size_%d_index_%d.keras' % (processed_sample, latent_dim, size, i)
      print(model_path)
      fout.write(model_path+"\n")
      collected_data = main(latent_dim ,model_path, processed_sample ,size, SAVE ) #latent_dim ,model_path, processed_sample ,size
      print(collected_data)
      fout.write(str(collected_data)+"\n"+"\n")
    fout.close()'''
    main()
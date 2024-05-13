# FROM : https://www.tensorflow.org/tutorials/generative/autoencoder 
#https://www.tensorflow.org/tutorials/keras/save_and_load 

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
import pickle
import os

work_dir = "/home/dina/Documents/PPI_WDLLM"

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
   
"""def mean_matrix_vec(interacting_prot, size): 
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
  return ca"""

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

def concatenate_pickle_ca(size, num_files):
    print ('Creating CA distance data input matrixes')
    # CA
    pickle_files_ca = [os.path.join(f"{work_dir}/Matrices_CA", file) for file in os.listdir(f"{work_dir}/Matrices_CA") if file.endswith('.pickle')]
        
    if num_files is not None:
        pickle_files_ca = pickle_files_ca[:num_files]
        
    concatenated_data_ca = []
    
    for file_name in pickle_files_ca:
        with open(file_name, 'rb') as g:
            data = np.array([pickle.load(g)])
            concatenated_data_ca.extend(data.reshape(-1, size, size))
    with open('concatenated_ca.pickle', 'wb') as g:
        pickle.dump(concatenated_data_ca, g)
    return concatenated_data_ca 

def concatenate_pickle_mean(size, num_files):
    # Mean
    print ('Creating Mean distance data input matrixes')
    pickle_files_mean = [os.path.join(f"{work_dir}/Matrices_Mean", file) for file in os.listdir(f"{work_dir}/Matrices_Mean") if file.endswith('.pickle')]
    if num_files is not None:
        pickle_files_mean = pickle_files_mean[:num_files]
    concatenated_data_mean = []
    for file_name in pickle_files_mean:
        with open(file_name, 'rb') as f:
            data_m = np.array([pickle.load(f)])
            concatenated_data_mean.extend(data_m.reshape(-1, size, size))
    with open('concatenated_mean.pickle', 'wb') as f:
        pickle.dump(concatenated_data_mean, f)

    return concatenated_data_mean

def main(latent_dim, model_name, processed_sample, size, SAVE, epoch):
  
  print('Running Autoencoder...')

  #feature_extraction.main(processed_sample, size)
  ranges= 2000
  int = 500
  loss_history = LossHistory()
  num_files = processed_sample
  me_a  = concatenate_pickle_mean(size, num_files)
  ca_a = concatenate_pickle_ca(size, num_files)
  
  dist_ca_train = np.array(me_a)
  dist_ca_test =  np.array(ca_a)
  dist_ca_train = dist_ca_train.astype('float32') / 255.
  dist_ca_test = dist_ca_test.astype('float32') / 255.
  print (f'Train input shape: {dist_ca_train.shape}')
  print (f'Test input shape: {dist_ca_test.shape}')
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
  print('Training ....')
  encoded_vectors_train = autoencoder.encoder(dist_ca_train).numpy()
  print(encoded_vectors_train.shape)
  #print(encoded_vectors_train)
  encoded_vectors_test = autoencoder.encoder(dist_ca_test).numpy()
  #print(dist_ca_train.shape)
  #print(encoded_vectors_train.shape)

  x, y , correlation, p_value, correlation2, p_value2  = spearman ( dist_ca_train, encoded_vectors_train , ranges , int )
  x_test, y_test , correlation_test, p_value_test, correlation2_test, p_value2_test  = spearman ( dist_ca_test, encoded_vectors_test , ranges , int )
  
  print("#################", model_name)
  print(correlation, p_value, correlation2, p_value2 )
  print(correlation_test, p_value_test, correlation2_test, p_value2_test )

  #plot(encoded_vectors_train, model_name)

  directory_path = '/home/dina/Documents/PPI_WDLLM/workdirautoencodermodelsbychain'
  model_name_with_path = os.path.join(directory_path, model_name)
  print('Saving model...')
  if SAVE:
    tf.saved_model.save(autoencoder, model_name_with_path)
  print('Done!')
  
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
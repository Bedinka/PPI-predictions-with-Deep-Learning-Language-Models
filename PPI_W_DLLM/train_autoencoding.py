# FROM : https://www.tensorflow.org/tutorials/generative/autoencoder 
#https://www.tensorflow.org/tutorials/keras/save_and_load 

import numpy as np
import pandas as pd
import tensorflow as tf
import os
import torch 
from scipy import stats
from random import randint
import pickle
import matplotlib.pyplot as plt

from autoencoder import Autoencoder
from sklearn.metrics import accuracy_score, precision_score, recall_score
from sklearn.model_selection import train_test_split
from tensorflow.keras import layers, losses, callbacks
from tensorflow.keras.models import Model
from keras.callbacks import Callback
from sklearn.decomposition import PCA



work_dir = "/home/dina/Documents/PPI_WDLLM"

    
class LossHistory(callbacks.Callback):
    def __init__(self):
        super().__init__()
        self.train_losses = []
        self.val_losses = []

    def on_epoch_end(self, epoch, logs=None):
        self.train_losses.append(logs.get('loss'))
        self.val_losses.append(logs.get('val_loss'))

def spearman(dist_ca_train,encoded_vectors_train, ranges , int  ):
  try:
    X = []
    Y = []
    random  = randint(1, 60)
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
    return X, Y , correlation, p_value, correlation2, p_value2 
  except:
     print('Couldnt perform Spearman')
  #print( correlation2, p_value2 )  
  

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

def split_data(data, train_size=0.8):
    total_size = len(data)
    train_end = int(train_size * total_size)
    
    indices = np.arange(total_size)
    np.random.shuffle(indices)
    
    train_idx = indices[:train_end].astype(int)
    test_idx = indices[train_end:].astype(int)
    
    train_data = data[train_idx]
    test_data = data[test_idx]
    
    return train_data, test_data

def concatenate_pickle(size, num_files, pickle_path):

    print('Create distance data input matrices as pickle')
    pickle_files = [os.path.join(pickle_path, file) for file in os.listdir(pickle_path) if file.endswith('.pickle')]
    
    if len(pickle_files) < num_files:
        print(f'Warning: Requested {num_files} files, but only {len(pickle_files)} available.')
        num_files = len(pickle_files)
    
    concatenated_data = []
    for file_name in pickle_files[:num_files]:
        with open(file_name, 'rb') as f:
            data_m = np.array([pickle.load(f)])
            concatenated_data.extend(np.array(data_m.reshape(-1, size, size)))
    return concatenated_data

def main(latent_dim, model_name, processed_sample, size, SAVE, epoch, batch_size , auto_dir , pickle_path):
  
  print('Running Autoencoder Training...')

  #TRAIN

  #model_name = 'dina_model_sample_%d_dim_%d_size_%d_epochs_%d.keras' % (processed_sample, latent_dim, size, epoch)
  ranges= 2000
  int = 500
  loss_history = LossHistory()
  num_files = processed_sample
  model_dir = os.path.join(work_dir, auto_dir)
  os.makedirs(model_dir, exist_ok=True)
  model_path = os.path.join(model_dir, model_name)
    
  data  = concatenate_pickle(size, num_files, pickle_path)
  train_data, test_data = split_data(data)
  #ca_a = concatenate_pickle_ca(size, num_files)
  
  dist_train = np.array(train_data)
  dist_test =  np.array(test_data)
  dist_train = dist_train.astype('float32') / 255.
  dist_test = dist_test.astype('float32') / 255.
  print (f'Train input shape: {dist_train.shape}')
  print (f'Test input shape: {dist_test.shape}')
  shape = dist_train.shape[1:]

  if os.path.exists(model_path):
    autoencoder = tf.keras.models.load_model(model_name, custom_objects={"Autoencoder": Autoencoder} )
    print("Loaded existing model.")
  else: 
    autoencoder = Autoencoder(latent_dim, shape)
    autoencoder.compile(optimizer='adam', loss=losses.MeanSquaredError())
    autoencoder.fit(dist_train, dist_train,
                      batch_size=batch_size,
                      epochs=epoch,
                      shuffle=True,
                      validation_data=(dist_test, dist_test),
                      callbacks=[loss_history])
  print('Training ....')
  encoded_vectors_train = autoencoder.encoder(dist_train).numpy()
  print(encoded_vectors_train.shape)
  #print(encoded_vectors_train)
  encoded_vectors_test = autoencoder.encoder(dist_test).numpy()
  #print(dist_ca_train.shape)
  #print(encoded_vectors_train.shape)

  x, y , correlation, p_value, correlation2, p_value2  = spearman ( dist_train, encoded_vectors_train , ranges , int )
  x_test, y_test , correlation_test, p_value_test, correlation2_test, p_value2_test  = spearman ( dist_ca_test, encoded_vectors_test , ranges , int )
  
  print("#################", model_name)
  print(correlation, p_value, correlation2, p_value2 )
  print(correlation_test, p_value_test, correlation2_test, p_value2_test )

  try:
      plot(encoded_vectors_train, model_name)
  except:
      print('Unable to plot.')

  model_name_with_path = os.path.join(auto_dir, model_name)
  print('Saving model...')
  if SAVE:
    tf.saved_model.save(autoencoder, model_name_with_path)

  print('Done!')

    
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
  return encoded_vectors_train, collected_data

if __name__ == "__main__":
  latent_dim = 2 
  processed_sample = 1600
  batch_size= 32
  size = 7
  SAVE = True
  epoch = 10
  auto_dir = './autoencoder_dina_models'
  model_name = 'dina_model_v2.keras'
  pickle_path = pickle_dir_ca = '/home/dina/Documents/PPI_WDLLM/Matrices_CA/train'

  main(latent_dim, model_name, processed_sample, size, SAVE, epoch, batch_size, auto_dir, pickle_path )

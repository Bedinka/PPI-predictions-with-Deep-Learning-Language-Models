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

    concatenated_data_ca = np.array(concatenated_data_ca)
    return  concatenated_data_ca

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

collected_data =[]

def main(model_name, size, unseen_data_path ):

    print('Loading trained model...')

    autoencoder = tf.keras.models.load_model(model_name, custom_objects={"Autoencoder": Autoencoder})
    if unseen_data_path is None:
        raise ValueError("Path to unseen data must be provided for testing.")
    with open(unseen_data_path, 'rb') as f:
        unseen_data = np.array([pickle.load(f)])
    unseen_data = unseen_data.reshape(-1, size, size).astype('float32') / 255.

    print('Testing with unseen data...')
    encoded_vectors_unseen = autoencoder.encoder(unseen_data).numpy()
    print(encoded_vectors_unseen.shape)
    
    ranges = 100
    interval = 10
    
    
    print("#################", model_name)
    try:
      x_unseen, y_unseen, correlation_unseen, p_value_unseen, correlation2_unseen, p_value2_unseen = spearman(unseen_data, encoded_vectors_unseen, ranges, interval)
      print(correlation_unseen, p_value_unseen, correlation2_unseen, p_value2_unseen)
      plot(encoded_vectors_unseen, model_name + "_unseen")
    except:
      print('Unable to Spearman or plot.')

    return encoded_vectors_unseen
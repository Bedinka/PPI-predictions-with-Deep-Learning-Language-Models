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

def split_data(file_list, train_size=0.8):
    print('Splitting Data into training and testing..')
    total_size = len(file_list)
    train_end = int(train_size * total_size)
    indices = np.arange(total_size)
    np.random.shuffle(indices)
    
    train_idx = indices[:train_end]
    test_idx = indices[train_end:]
    
    train_files = [file_list[i] for i in train_idx]
    test_files = [file_list[i] for i in test_idx]

    print(f'Training size: {len(train_files):,}')
    print(f'Test size: {len(test_files):,}')

    return train_files, test_files

def concatenate_pickle(size, num_files, pickle_path):

    print('Create distance data input matrices as pickle')
    pickle_files = [os.path.join(pickle_path, file) for file in os.listdir(pickle_path) if file.endswith('.pickle')]
    concatenated_pickle = 'concatenated_subs.pickle'
    if num_files is not None:
          pickle_files = pickle_files[:num_files]
    concatenated_data = []
    for file_name in pickle_files:
        with open(file_name, 'rb') as f:
            data_m = np.array([pickle.load(f)])
            concatenated_data.extend(data_m.reshape(-1, size, size))
    with open(concatenated_pickle, 'wb') as f:
      pickle.dump(concatenated_data, f)
    concatenated_data = np.array(concatenated_data)
    return  concatenated_data


def create_tf_dataset(pickle_files, batch_size, size):
    print('Loading in the pickle files..')
    def data_generator():
        for file_name in pickle_files:
            with open(file_name, 'rb') as f:
                data = np.array(pickle.load(f)).reshape(-1, size, size).astype('float32') / 255.
                num_batches = len(data) // batch_size
                for i in range(num_batches):
                    batch_data = data[i * batch_size:(i + 1) * batch_size]
                    yield batch_data, batch_data

    dataset = tf.data.Dataset.from_generator(
        data_generator,
        output_signature=(
            tf.TensorSpec(shape=(batch_size, size, size), dtype=tf.float32),
            tf.TensorSpec(shape=(batch_size, size, size), dtype=tf.float32)
        )
    )
    return dataset.prefetch(tf.data.experimental.AUTOTUNE)

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'


def main(latent_dim, model_name, processed_sample, size, SAVE, epoch, batch_size , auto_dir , pickle_path):
  
  print('Running Autoencoder Training...')

  #TRAIN

  #model_name = 'dina_model_sample_%d_dim_%d_size_%d_epochs_%d.keras' % (processed_sample, latent_dim, size, epoch)
  ranges= 2000
  int = 500
  loss_history = LossHistory()
  num_files = processed_sample
  
  pickle_files = [os.path.join(pickle_path, file) for file in os.listdir(pickle_path) if file.endswith('.pickle')][:processed_sample]
  load_batch_size = 128
  train_files, test_files = split_data(pickle_files)

  train_dataset = create_tf_dataset(train_files, load_batch_size, size)
  test_dataset = create_tf_dataset(test_files, load_batch_size, size)
  os.makedirs(auto_dir, exist_ok=True)
  model_path = os.path.join(auto_dir, model_name)

  if os.path.exists(model_path):
      autoencoder = tf.keras.models.load_model(model_path, custom_objects={"Autoencoder": Autoencoder})
      print("Loaded existing model.")
  else:
      autoencoder = Autoencoder(latent_dim, (size, size))
      autoencoder.compile(optimizer='adam', loss=losses.MeanSquaredError())
      autoencoder.fit(train_dataset,
                      epochs=epoch,
                      validation_data=test_dataset,
                      callbacks=[loss_history],
                      batch_size=batch_size, 
                      verbose=0)
 
  print('Training ....')
  train_data = np.concatenate([x for x, _ in train_dataset], axis=0)
  test_data = np.concatenate([x for x, _ in test_dataset], axis=0)
  
  encoded_vectors_train = autoencoder.encoder(train_data).numpy()
  print(encoded_vectors_train.shape)
  #print(encoded_vectors_train)
  encoded_vectors_test = autoencoder.encoder(test_data).numpy()
  #print(dist_ca_train.shape)
  #print(encoded_vectors_train.shape)

  x, y , correlation, p_value, correlation2, p_value2  = spearman ( train_data, encoded_vectors_train , ranges , int )
  x_test, y_test , correlation_test, p_value_test, correlation2_test, p_value2_test  = spearman ( test_data, encoded_vectors_test , ranges , int )
  
  print("#################", model_name)
  print(correlation, p_value, correlation2, p_value2 )
  print(correlation_test, p_value_test, correlation2_test, p_value2_test )

  try:
      plot(encoded_vectors_train, model_name)
  except:
      print('Unable to plot.')

  model_name_with_path = os.path.join(auto_dir, model_name)
  
  if SAVE:
    print('Saving model...')
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
  processed_sample = 30
  batch_size= 32
  size = 7
  SAVE = True
  epoch = 10
  auto_dir = '/home/dina/Documents/PPI_WDLLM/autoencoder_dina_models'
  model_name = 'dina_model_v2.keras'
  pickle_path =  '/home/dina/Documents/PPI_WDLLM/Matrices_CA/train'

  main(latent_dim, model_name, processed_sample, size, SAVE, epoch, batch_size, auto_dir, pickle_path )

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
    dist_matrix = np.array([protein.residues, [protein.distance_matrices_CA,[ protein.distance_matrices_mean]]], dtype= object)
    print(dist_matrix.shape)
    tensor = tf.constant(dist_matrix, dtype = object , name = 'tensor1')
    (dist_ca_train, _), (dist_ca_test, _) = dist_matrix

    """
    Traceback (most recent call last):
  File "/home/pc550/Documents/PPI_W_DLLM/YangLabIntern/PPI_W_DLLM/data_conversion.py", line 16, in <module>
    dist_matrix = np.array([protein.residues, [protein.distance_matrices_CA,[ protein.distance_matrices_mean]]])
ValueError: setting an array element with a sequence. The requested array has an inhomogeneous shape after 1 dimensions. The detected shape was (2,) + inhomogeneous part.

Traceback (most recent call last):
  File "/home/pc550/Documents/PPI_W_DLLM/YangLabIntern/PPI_W_DLLM/data_conversion.py", line 18, in <module>
    tensor = tf.convert_to_tensor(dist_matrix, dtype=object, name='tensor1')
  File "/home/pc550/miniconda3/envs/train/lib/python3.8/site-packages/tensorflow/python/util/traceback_utils.py", line 153, in error_handler
    raise e.with_traceback(filtered_tb) from None
  File "/home/pc550/miniconda3/envs/train/lib/python3.8/site-packages/tensorflow/python/framework/constant_op.py", line 98, in convert_to_eager_tensor
    return ops.EagerTensor(value, ctx.device_name, dtype)
ValueError: Failed to convert a NumPy array to a Tensor (Unsupported object type list).


    (dist_ca_ab_train, _)(dist_ca_ab_test, _)= protein.distance_matrices_CA_AB
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
latent_dim = 128
autoencoder = Autoencoder(latent_dim, shape)

autoencoder.compile(optimizer='adam', loss=losses.MeanSquaredError())

autoencoder.fit(dist_ca_train, dist_ca_train,
                epochs=10,
                shuffle=True,
                validation_data=(dist_ca_test, dist_ca_test))


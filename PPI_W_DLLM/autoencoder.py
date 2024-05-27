import tensorflow as tf
from tensorflow.keras.models import Model

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
    
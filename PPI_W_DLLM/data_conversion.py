import torch
import numpy as np
import feature_extraction
from tensorflow.keras.layers import AveragePooling2D
from tensorflow.keras.layers import Input, Concatenate, Conv2D, Flatten, Dense, Reshape, GlobalAveragePooling2D
from tensorflow.keras.models import Model

interacting_prot, chains_CA = feature_extraction.main()

def create_autoencoder_with_pooling(input_shape, num_classes, latent_dim):
    matrix_input = Input(shape=input_shape)

    class_input = Input(shape=(num_classes,))

    pooled_input = AveragePooling2D(pool_size=(2, 2))(matrix_input)

    concatenated_input = Concatenate()([pooled_input, class_input])

    encoder_output = GlobalAveragePooling2D()(Conv2D(32, (3, 3), activation='relu')(concatenated_input))

    decoder_input = Input(shape=(latent_dim,))
    decoder_concatenated_input = Concatenate()([decoder_input, class_input])
    decoder_output = Dense(np.prod(input_shape), activation='sigmoid')(decoder_concatenated_input)
    decoded_output = Reshape(input_shape)(decoder_output)


    encoder = Model(inputs=[matrix_input, class_input], outputs=encoder_output)
    decoder = Model(inputs=[decoder_input, class_input], outputs=decoded_output)
    autoencoder = Model(inputs=[matrix_input, class_input], outputs=decoder(encoder([matrix_input, class_input])))

    return autoencoder

def extract_matrices_and_classes(chains_CA):
    all_matrices = []
    all_classes = []
    for chain_id, chain in chains_CA.items():
        all_matrices.extend(chain.distance_matrices_CA_AB)
        all_matrices.extend(chain.distance_matrices_mean_AB)
        all_matrices.extend(chain.submatrices)
        all_classes.extend(chain.__dict__.values())
    return all_matrices, all_classes
 
def main ():
    all_matrices, all_classes = extract_matrices_and_classes(chains_CA)
    max_height = max_width = max_channels = 0
    for matrix in all_matrices:
            height, width = matrix.row, matrix.col
            max_height = max(max_height, height)
            max_width = max(max_width, width)
            

    # input shape using maximum shape
    input_shape = (max_height, max_width, max_channels)
    print(input_shape)
    # number of classes and latent dimension for autoencoder
    num_classes = len(interacting_prot)         
    latent_dim = num_classes     
    print(num_classes)                      
    # autoencoder model with pooling
    autoencoder = create_autoencoder_with_pooling(input_shape, num_classes, latent_dim)

    autoencoder.fit([all_matrices, all_classes], all_matrices, epochs=epochs, batch_size=batch_size)

    # fixed-size vector representation from the encoder for each PDB file
    encoded_vectors = autoencoder.get_layer("encoder").predict([all_matrices, all_classes])

if __name__ == "__main__":
    main()
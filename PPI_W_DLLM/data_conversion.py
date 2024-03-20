import numpy as np
import feature_extraction
from tensorflow.keras.layers import AveragePooling2D, Masking, LSTM, Concatenate, Dense, Reshape
from tensorflow.keras.layers import Input, Concatenate, Conv2D, Flatten, Dense, Reshape, GlobalAveragePooling2D
from tensorflow.keras.models import Model

interacting_prot, chains_CA = feature_extraction.main()


def create_autoencoder_with_pooling(input_shape, num_classes, latent_dim):
    class_input = Input(shape=(num_classes,))
    chain1_input = Input(shape=(input_shape[0]))  
    chain2_input = Input(shape=(input_shape[0])) 

    pooled_input_chain1 = AveragePooling2D(pool_size=(2, 2))(chain1_input)
    pooled_input_chain2 = AveragePooling2D(pool_size=(2, 2))(chain2_input)

    concatenated_input = Concatenate()([pooled_input_chain1, pooled_input_chain2, class_input])

    encoder_output = GlobalAveragePooling2D()(Conv2D(32, (3, 3), activation='relu')(concatenated_input))

    decoder_input = Input(shape=(latent_dim,))
    decoder_concatenated_input = Concatenate()([decoder_input, class_input])
    decoder_output = Dense(np.prod(input_shape), activation='sigmoid')(decoder_concatenated_input)
    decoded_output = Reshape(input_shape)(decoder_output)

    encoder = Model(inputs=[chain1_input, chain2_input, class_input], outputs=encoder_output)
    decoder = Model(inputs=[decoder_input, class_input], outputs=decoded_output)
    autoencoder = Model(inputs=[chain1_input, chain2_input, class_input], outputs=decoder(encoder([chain1_input, chain2_input, class_input])))

    return autoencoder

def extract_classes(chains_CA):
    all_classes = []
    for chain_id, chain in chains_CA.items():
        all_classes.append(chain.__dict__.values())
    return all_classes

def main():
    all_classes = extract_classes(chains_CA)
    # Determine input shapes dynamically

    input_shape =  100
    max_channels = 17
    input_shapes = (input_shape, max_channels)
    num_classes = 2
    latent_dim = 2

    # Create autoencoder model
    autoencoder = create_autoencoder_with_pooling(input_shapes, num_classes, latent_dim)
    autoencoder.compile(optimizer='adam', loss='mse')
    autoencoder.fit(interacting_prot)

    # Extract encoded vectors
    encoded_vectors = autoencoder.layers[1].predict()

if __name__ == "__main__":
    main()

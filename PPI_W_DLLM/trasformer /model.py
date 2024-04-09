import torch 
import torch.nn as nn
import math

class LayerNormalization(nn.Module):
    # item in the sequence, we calculate for each item independently the mean and variance , calculate ne wvalues with theit own mean and variance 
    def __init__(self, features: int, eps:float=10**-6) -> None:
        super().__init__()
        self.eps = eps # epsilon needed in the normalization formula , so the x value doesnt get too small or  big , if sigma changes, numeric stability 
        self.alpha = nn.Parameter(torch.ones(features)) # alpha is a learnable parameter, multiplied ( nn makes the parameter learnable)
        self.bias = nn.Parameter(torch.zeros(features)) # bias is a learnable parameter, additive

    def forward(self, x):
        # x: (batch, seq_len, hidden_size)
         # Keep the dimension for broadcasting, everything after the bacth , keepdim : mean canceles the dimention which is it applied but we want to keep it 
        mean = x.mean(dim = -1, keepdim = True) # (batch, seq_len, 1)
        # Keep the dimension for broadcasting, 
        std = x.std(dim = -1, keepdim = True) # (batch, seq_len, 1)
        # eps is to prevent dividing by zero or when std is very small, by the FROMULA
        return self.alpha * (x - mean) / (std + self.eps) + self.bias

class FeedForwardBlock(nn.Module):
    # fully connected layer that the model uses both encoder and decoder parts, basically two matrixes multiplied by x with relu inbetween with a bias
    def __init__(self, d_model: int, d_ff: int, dropout: float) -> None: # d_model = 512 , d_ff = 2048 ???? WHY 
        super().__init__()
        self.linear_1 = nn.Linear(d_model, d_ff) # w1 and b1 , first matrix 
        self.dropout = nn.Dropout(dropout)
        self.linear_2 = nn.Linear(d_ff, d_model) # w2 and b2, second matrix 

    def forward(self, x):
        # (batch, seq_len, d_model) --> (batch, seq_len, d_ff) --> (batch, seq_len, d_model)
        return self.linear_2(self.dropout(torch.relu(self.linear_1(x))))
    

#Embedding 
class InputEmbedding(nn.Module):

    def __init__(self, d_model: int, vocab_size: int):
        super().__init__()
        self.d_model= d_model
        self.vocab_size = vocab_size
        # Embending is usually taken a value and create a FIXED SIZED (d_modell size ) vector, mapping between numbers and the vector of size
        self.embedding = nn.Embedding( (vocab_size, d_model))
    3
    def forward(self, x):
        return self.embedding(x) * math.sqrt(self.d_model) # taking the embedding,  multiply it by squeare root of d_modell , Pytorch own layer which does the mapping for us to create 512 

# Positional encoding : mapped by a list of vectors by the embedding 
# adding a vector of the same size as the embedding to th embedding vector , which containes the position of item in the sequence 

class PositionalEncoding (nn.Module):
    def __init(self, d_model: int, seq_len: int , dropout :  float):  #one vector for each postion 
        super.__init__()
        self.d_modell= d_model  #modell size defines the size of the vector 
        self.seq_len = seq_len   #max sequence length cause wee need to give one vector for each position 
        self.dropout = nn.Dropout(dropout) # dropout make the modell less overfitt

        #matrix of shape seqlenght to d_modell, vectors of d_modell size but seq number of them 
        pe= torch.zeros(seq_len, d_model) # creating matrix of shape (seq, d_model)
        position = torch.arrange(0, seq_len, dtype = torch.float).unsqueeze(1) #creating vector of shape seq_len , (seq_len,1) 
        
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model)) # diviator of the formula PE , (d_model / 2)
        # Apply sine to even indices
        pe[:, 0::2] = torch.sin(position * div_term) # sin(position * (10000 ** (2i / d_model))
        # Apply cosine to odd indices
        pe[:, 1::2] = torch.cos(position * div_term) # cos(position * (10000 ** (2i / d_model))
        # Add a batch dimension to the positional encoding, so it can be aplied to the whole sequence , now is (seq, d_model), ad a new dim 
        pe = pe.unsqueeze(0) # (1, seq_len, d_model), ad a new dimensions
        # Register the positional encoding as a buffer
        self.register_buffer('pe', pe) # buffer keep this values saved when the modell is also saved , should be registered as a buffer 

    def forward(self, x): # add this positional encoding to each item if the sequence 
        x = x + (self.pe[:, :x.shape[1], :]).requires_grad_(False) # (batch, seq_len, d_model) , they will be always fixed and not learned(requires_grad)
        return self.dropout(x)

class ResidualConnection(nn.Module):

        def __init__(self, features: int, dropout: float) -> None: 
            super().__init__()
            self.dropout = nn.Dropout(dropout)
            self.norm = LayerNormalization(features)
    
        def forward(self, x, sublayer):
            return x + self.dropout(sublayer(self.norm(x)))

class MultiHeadAttentionBlock(nn.Module):
    # takes the input of the encoder and uses it 3 times, Query, Keys, Values , in the encoder these matrixes will all be the same , not in the decoder tho !
    # (Q, K , V ) X ( Wq, Wk, Wv) = (Q', K', V') --> H,  smaller matrixes for each head , splitting these by the embedding dimension, each head will have acces to the whole sequence , and attention to aplied to all of them, then concatenate them all on the h dimension and multiplie by Wo = MH-A(seq, d_model)
    # there is a batch dimesion dnt forget 
    def __init__(self, d_model: int, h: int, dropout: float) -> None:# h : number of heads
        super().__init__()
        self.d_model = d_model # Embedding vector size
        self.h = h # Number of heads
        # Make sure d_model is divisible by h
        assert d_model % h == 0, "d_model is not divisible by h, dividing embedding vector into H heads , which mneed to divide equally the same domain of the same vector representing the mbedding into equal matrices to ech "

        self.d_k = d_model // h # Dimension of vector seen by each head
        self.w_q = nn.Linear(d_model, d_model, bias=False) # Wq ( d_m, d_m)! this will give back after the multiplication with (sed, d_m) = Q' (seq, d_m)
        self.w_k = nn.Linear(d_model, d_model, bias=False) # Wk
        self.w_v = nn.Linear(d_model, d_model, bias=False) # Wv
        self.w_o = nn.Linear(d_model, d_model, bias=False) # Wo ( h* dv = d_m  , d_m)
        self.dropout = nn.Dropout(dropout)

    @staticmethod # i can call this function without having an instance of this class 
    def attention(query, key, value, mask, dropout: nn.Dropout): # mask is when we dont want items to interract with each other we mask them  , in our case the self interaction are 0
        #D-k last dim of the Q, K, V 
        d_k = query.shape[-1]
        # Just apply the formula from the paper
        # (batch, h, seq_len, d_k) --> (batch, h, seq_len, seq_len)
        attention_scores = (query @ key.transpose(-2, -1)) / math.sqrt(d_k) # @ add sign is matrix multiplication in pytorch , transpose the last two dimensions se by d_k and D_k by seq_len 
        if mask is not None:
            # Write a very low value (indicating -inf) to the positions where mask == 0
            attention_scores.masked_fill_(mask == 0, -1e9)
        attention_scores = attention_scores.softmax(dim=-1) # (batch, h, seq_len, seq_len) # Apply softmax
        if dropout is not None:
            attention_scores = dropout(attention_scores)
        # (batch, h, seq_len, seq_len) --> (batch, h, seq_len, d_k)
        # return attention scores which can be used for visualization
        return (attention_scores @ value), attention_scores

    def forward(self, q, k, v, mask): # during the softmax we create a sequence by sequence matrix , we can mask values in this matrix before aplying the softmax 
        query = self.w_q(q) # (batch, seq_len, d_model) --> (batch, seq_len, d_model) Q'
        key = self.w_k(k) # (batch, seq_len, d_model) --> (batch, seq_len, d_model) K'
        value = self.w_v(v) # (batch, seq_len, d_model) --> (batch, seq_len, d_model)V'

        # (batch, seq_len, d_model) --> (batch, seq_len, h, d_k) --> (batch, h, seq_len, d_k)  we want each head to watch (seq-len, d_k) , so full sequence  whit a small part of embedding
        query = query.view(query.shape[0], query.shape[1], self.h, self.d_k).transpose(1, 2) # .view in pytorch method , which keeps the batch dimension ( we want to split the  embedding not the sequence into H parts) , in the 2. dim we want to keep the sequence, d_model split into two smaller dim which h by d_k , gives you back the d_m shape 
        # we transpose cause H dimension is prefered, so each item can see the whole sequence
        key = key.view(key.shape[0], key.shape[1], self.h, self.d_k).transpose(1, 2) # 
        value = value.view(value.shape[0], value.shape[1], self.h, self.d_k).transpose(1, 2)

        # Calculate attention
        x, self.attention_scores = MultiHeadAttentionBlock.attention(query, key, value, mask, self.dropout) # we need the output and the attention score ( out of softmax)
        
        # Combine all the heads together
        # (batch, h, seq_len, d_k) --> (batch, seq_len, h, d_k) --> (batch, seq_len, d_model)
        x = x.transpose(1, 2).contiguous().view(x.shape[0], -1, self.h * self.d_k)

        # Multiply by Wo
        # (batch, seq_len, d_model) --> (batch, seq_len, d_model)  
        return self.w_o(x)

class EncoderBlock(nn.Module):


    def __init__(self, features: int, self_attention_block: MultiHeadAttentionBlock, feed_forward_block: FeedForwardBlock, dropout: float) -> None:
        super().__init__()
        self.self_attention_block = self_attention_block
        self.feed_forward_block = feed_forward_block
        self.residual_connections = nn.ModuleList([ResidualConnection(features, dropout) for _ in range(2)])

    def forward(self, x, src_mask):
        x = self.residual_connections[0](x, lambda x: self.self_attention_block(x, x, x, src_mask))
        x = self.residual_connections[1](x, self.feed_forward_block)
        return x
    
class Encoder(nn.Module):

    def __init__(self, features: int, layers: nn.ModuleList) -> None:
        super().__init__()
        self.layers = layers
        self.norm = LayerNormalization(features)

    def forward(self, x, mask):
        for layer in self.layers:
            x = layer(x, mask)
        return self.norm(x)

class DecoderBlock(nn.Module):

    def __init__(self, features: int, self_attention_block: MultiHeadAttentionBlock, cross_attention_block: MultiHeadAttentionBlock, feed_forward_block: FeedForwardBlock, dropout: float) -> None:
        super().__init__()
        self.self_attention_block = self_attention_block
        self.cross_attention_block = cross_attention_block
        self.feed_forward_block = feed_forward_block
        self.residual_connections = nn.ModuleList([ResidualConnection(features, dropout) for _ in range(3)])

    def forward(self, x, encoder_output, src_mask, tgt_mask):
        x = self.residual_connections[0](x, lambda x: self.self_attention_block(x, x, x, tgt_mask))
        x = self.residual_connections[1](x, lambda x: self.cross_attention_block(x, encoder_output, encoder_output, src_mask))
        x = self.residual_connections[2](x, self.feed_forward_block)
        return x
    
class Decoder(nn.Module):

    def __init__(self, features: int, layers: nn.ModuleList) -> None:
        super().__init__()
        self.layers = layers
        self.norm = LayerNormalization(features)

    def forward(self, x, encoder_output, src_mask, tgt_mask):
        for layer in self.layers:
            x = layer(x, encoder_output, src_mask, tgt_mask)
        return self.norm(x)

class ProjectionLayer(nn.Module):

    def __init__(self, d_model, vocab_size) -> None:
        super().__init__()
        self.proj = nn.Linear(d_model, vocab_size)

    def forward(self, x) -> None:
        # (batch, seq_len, d_model) --> (batch, seq_len, vocab_size)
        return self.proj(x)
    
class Transformer(nn.Module):

    def __init__(self, encoder: Encoder, decoder: Decoder, src_embed: InputEmbeddings, tgt_embed: InputEmbeddings, src_pos: PositionalEncoding, tgt_pos: PositionalEncoding, projection_layer: ProjectionLayer) -> None:
        super().__init__()
        self.encoder = encoder
        self.decoder = decoder
        self.src_embed = src_embed
        self.tgt_embed = tgt_embed
        self.src_pos = src_pos
        self.tgt_pos = tgt_pos
        self.projection_layer = projection_layer

    def encode(self, src, src_mask):
        # (batch, seq_len, d_model)
        src = self.src_embed(src)
        src = self.src_pos(src)
        return self.encoder(src, src_mask)
    
    def decode(self, encoder_output: torch.Tensor, src_mask: torch.Tensor, tgt: torch.Tensor, tgt_mask: torch.Tensor):
        # (batch, seq_len, d_model)
        tgt = self.tgt_embed(tgt)
        tgt = self.tgt_pos(tgt)
        return self.decoder(tgt, encoder_output, src_mask, tgt_mask)
    
    def project(self, x):
        # (batch, seq_len, vocab_size)
        return self.projection_layer(x)
    
def build_transformer(src_vocab_size: int, tgt_vocab_size: int, src_seq_len: int, tgt_seq_len: int, d_model: int=512, N: int=6, h: int=8, dropout: float=0.1, d_ff: int=2048) -> Transformer:
    # Create the embedding layers
    src_embed = InputEmbeddings(d_model, src_vocab_size)
    tgt_embed = InputEmbeddings(d_model, tgt_vocab_size)

    # Create the positional encoding layers
    src_pos = PositionalEncoding(d_model, src_seq_len, dropout)
    tgt_pos = PositionalEncoding(d_model, tgt_seq_len, dropout)
    
    # Create the encoder blocks
    encoder_blocks = []
    for _ in range(N):
        encoder_self_attention_block = MultiHeadAttentionBlock(d_model, h, dropout)
        feed_forward_block = FeedForwardBlock(d_model, d_ff, dropout)
        encoder_block = EncoderBlock(d_model, encoder_self_attention_block, feed_forward_block, dropout)
        encoder_blocks.append(encoder_block)

    # Create the decoder blocks
    decoder_blocks = []
    for _ in range(N):
        decoder_self_attention_block = MultiHeadAttentionBlock(d_model, h, dropout)
        decoder_cross_attention_block = MultiHeadAttentionBlock(d_model, h, dropout)
        feed_forward_block = FeedForwardBlock(d_model, d_ff, dropout)
        decoder_block = DecoderBlock(d_model, decoder_self_attention_block, decoder_cross_attention_block, feed_forward_block, dropout)
        decoder_blocks.append(decoder_block)
    
    # Create the encoder and decoder
    encoder = Encoder(d_model, nn.ModuleList(encoder_blocks))
    decoder = Decoder(d_model, nn.ModuleList(decoder_blocks))
    
    # Create the projection layer
    projection_layer = ProjectionLayer(d_model, tgt_vocab_size)
    
    # Create the transformer
    transformer = Transformer(encoder, decoder, src_embed, tgt_embed, src_pos, tgt_pos, projection_layer)
    
    # Initialize the parameters
    for p in transformer.parameters():
        if p.dim() > 1:
            nn.init.xavier_uniform_(p)
    
    return transformer
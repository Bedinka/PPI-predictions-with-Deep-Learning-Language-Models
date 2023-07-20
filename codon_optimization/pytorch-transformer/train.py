import re
from model import build_transformer
from dataset import BilingualDataset, causal_mask
from datasets import Dataset
    
from config import get_config, get_weights_file_path

import torchtext.datasets as datasets
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, random_split
from torch.optim.lr_scheduler import LambdaLR

import warnings
from tqdm import tqdm
import os
from pathlib import Path

# Huggingface datasets and tokenizers
from datasets import load_dataset
from tokenizers import Tokenizer
from tokenizers.models import WordLevel
from tokenizers.trainers import WordLevelTrainer
from tokenizers.pre_tokenizers import Whitespace

import torchmetrics
from torch.utils.tensorboard import SummaryWriter

codon_dic = {"GCT":"A","GCC":"A","GCA":"A","GCG":"A","TGT":"C","TGC":"C","GAA":"E","GAG":"E","GAT":"D","GAC":"D","GGT":"G","GGC":"G","GGA":"G","GGG":"G","TTT":"F","TTC":"F","ATT":"I","ATC":"I","ATA":"I","CAT":"H","CAC":"H","AAA":"K","AAG":"K","ATG":"M","CTT":"L","CTC":"L","CTA":"L","CTG":"L","TTA":"L","TTG":"L","AAT":"N","AAC":"N","CAA":"Q","CAG":"Q","CCT":"P","CCC":"P","CCA":"P","CCG":"P","TCT":"S","TCC":"S","TCA":"S","TCG":"S","AGT":"S","AGC":"S","CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R","ACT":"T","ACC":"T","ACA":"T","ACG":"T","TGG":"W","GTT":"V","GTC":"V","GTA":"V","GTG":"V","TAT":"Y","TAC":"Y"}

aa_dic = { "I" : ["ATT", "ATC", "ATA"], "L" : ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"], "V" : ["GTT", "GTC", "GTA", "GTG"], "F" : ["TTT", "TTC"], "M" : ["ATG"], "C" : ["TGT", "TGC"], "A" : ["GCT", "GCC", "GCA", "GCG"], "G" : ["GGT", "GGC", "GGA", "GGG"], "P" : ["CCT", "CCC", "CCA", "CCG"], "T" : ["ACT", "ACC", "ACA", "ACG"], "S" : ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], "Y" : ["TAT", "TAC"], "W" : ["TGG"], "Q" : ["CAA", "CAG"], "N" : ["AAT", "AAC"], "H" : ["CAT", "CAC"], "E" : ["GAA", "GAG"], "D" : ["GAT", "GAC"], "K" : ["AAA", "AAG"], "R" : ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"] }


def greedy_decode(model, source, source_mask, tokenizer_src, tokenizer_tgt, max_len, device):
    sos_idx = tokenizer_tgt.token_to_id('[SOS]')
    eos_idx = tokenizer_tgt.token_to_id('[EOS]')

    # Precompute the encoder output and reuse it for every step
    encoder_output = model.encode(source, source_mask)
    # Initialize the decoder input with the sos token
    decoder_input = torch.empty(1, 1).fill_(sos_idx).type_as(source).to(device)
    while True:
        if decoder_input.size(1) == max_len:
            break

        # build mask for target
        decoder_mask = causal_mask(decoder_input.size(1)).type_as(source_mask).to(device)

        # calculate output
        out = model.decode(encoder_output, source_mask, decoder_input, decoder_mask)

        # get next token
        prob = model.project(out[:, -1])
        _, next_word = torch.max(prob, dim=1)
        decoder_input = torch.cat(
            [decoder_input, torch.empty(1, 1).type_as(source).fill_(next_word.item()).to(device)], dim=1
        )

        if next_word == eos_idx:
            break

    return decoder_input.squeeze(0)

def run_validation(model, validation_ds, tokenizer_src, tokenizer_tgt, max_len, device, print_msg, global_step, writer, num_examples=2):
    model.eval()
    count = 0

    source_texts = []
    expected = []
    predicted = []

    try:
        # get the console window width
        with os.popen('stty size', 'r') as console:
            _, console_width = console.read().split()
            console_width = int(console_width)
    except:
        # If we can't get the console width, use 80 as default
        console_width = 80

    with torch.no_grad():
        for batch in validation_ds:
            count += 1
            encoder_input = batch["encoder_input"].to(device) # (b, seq_len)
            encoder_mask = batch["encoder_mask"].to(device) # (b, 1, 1, seq_len)

            # check that the batch size is 1
            assert encoder_input.size(
                0) == 1, "Batch size must be 1 for validation"

            model_out = greedy_decode(model, encoder_input, encoder_mask, tokenizer_src, tokenizer_tgt, max_len, device)

            source_text = batch["src_text"][0]
            target_text = batch["tgt_text"][0]
            model_out_text = tokenizer_tgt.decode(model_out.detach().cpu().numpy())
            print_msg(f"pred_length: {len(model_out_text)}; target_length:{len(target_text)}")
            model_out_text += " [EOS]"

            source_texts.append(source_text)
            expected.append(target_text)
            predicted.append(model_out_text)
            
            # Print the source, target and model output
            print_msg('-'*console_width)
            print_msg(f"{f'SOURCE: ':>12}{source_text}")
            print_msg(f"{f'TARGET: ':>12}{target_text}")
            print_msg(f"{f'PREDICTED: ':>12}{model_out_text}")

            mRNA_original = target_text.replace(" ", "")
            mRNA_predicted = model_out_text.replace(" ", "")
            aa_original = translate( mRNA_original )
            aa_predicted = translate( mRNA_predicted )    
            
            print_msg(f"{f'TARGET AAs: ':>12}{aa_original}")
            print_msg(f"{f'PREDICTED AAs: ':>12}{aa_predicted}")
            print_msg(f"{f'SIMILARITY SCORE: ':>12}{similarity( aa_original, aa_predicted )}")

            import os

            filename = "sim_scores.txt"
            counter = 1

            while os.path.exists(filename):
                name, ext = os.path.splitext(filename)
                
                new_filename = f"{name}_{counter}{ext}"
                counter += 1
                
                filename = new_filename

            with open(filename, "w") as f:
                f.write(str(similarity( aa_original, aa_predicted )) + '\n')


            if count == num_examples:
                print_msg('-'*console_width)
                break
        
    
    if writer:
        # Evaluate the character error rate
        # Compute the char error rate 
        metric = torchmetrics.CharErrorRate()
        cer = metric(predicted, expected)
        writer.add_scalar('validation cer', cer, global_step)
        writer.flush()

        # Compute the word error rate
        metric = torchmetrics.WordErrorRate()
        wer = metric(predicted, expected)
        writer.add_scalar('validation wer', wer, global_step)
        writer.flush()

        # Compute the BLEU metric
        metric = torchmetrics.BLEUScore()
        bleu = metric(predicted, expected)
        writer.add_scalar('validation BLEU', bleu, global_step)
        writer.flush()

def get_all_sentences(ds, lang):
    for item in ds:
        yield item['translation'][lang]

def get_or_build_tokenizer(config, ds, lang):
    tokenizer_path = Path(config['tokenizer_file'].format(lang))
    if not Path.exists(tokenizer_path):
        # Most code taken from: https://huggingface.co/docs/tokenizers/quicktour
        tokenizer = Tokenizer(WordLevel(unk_token="[UNK]"))
        tokenizer.pre_tokenizer = Whitespace()
        trainer = WordLevelTrainer(special_tokens=["[UNK]", "[PAD]", "[SOS]", "[EOS]"], min_frequency=2)
        tokenizer.train_from_iterator(get_all_sentences(ds, lang), trainer=trainer)
        tokenizer.save(str(tokenizer_path))
    else:
        tokenizer = Tokenizer.from_file(str(tokenizer_path))
    return tokenizer

def build_tokenizer_from_data(data, tokenizer_path):
    # data = ["AAAAGASDAS", "GGGTHA", ...] # array of string
    # Most code taken from: https://huggingface.co/docs/tokenizers/quicktour
    # https://huggingface.co/docs/tokenizers/python/v0.10.0/tutorials/python/training_from_memory.html
    tokenizer = Tokenizer(WordLevel(unk_token="[UNK]"))
    tokenizer.pre_tokenizer = Whitespace()
    trainer = WordLevelTrainer(special_tokens=["[UNK]", "[PAD]", "[SOS]", "[EOS]"], min_frequency=2)
    tokenizer.train_from_iterator(data, trainer=trainer)
    tokenizer.save(str(tokenizer_path))

MAX_SEQ_LEN = 100 # temporary

def trans_inititation_split(seq):
    m = re.search('(.*)ATG(.{39})', seq, re.IGNORECASE)

    if m.group(1) == None:
        pre = "N" * 39

    while m.group(1) != None:
        pre = m.group(1)
        while len(pre) < 39:
            pre = "N" + m.group(1)
        while len(pre) > 39:
            pre = pre[len(pre)-39:]

    trans_init_seq = "ATG".join(pre, m.group(2))

    rest = seq[m.pos+39:]
    return trans_init_seq, rest

def trans_termination_split(seq):
    m = re.search('(.{39})[TGA,TGG,TAA](.{39})', seq, re.IGNORECASE)
    trans_term_seq = "TAA".join(m.group(1,2))
    rest = seq[:m.pos-39]
    return trans_term_seq, rest

def middle_split(seq):
    middlepos = len(seq)//2
    middle = seq[middlepos-39:middlepos+39]
    return middle

def aa_start(aa):
    return aa[:13]

def split_word_amino_acid_seq(seq):
    return " ".join(seq)

def split_word_mRNA(mRNA):
    return " ".join( [ mRNA[i:i+3] for i in range(0, len(mRNA), 3) ] )

def load_amino_acid_data(input_filepath):
    data = []
    f = open( input_filepath )
    for line in f.readlines():
        if line[0] == ">": continue
        seq = line.strip()
        data.append( split_word_amino_acid_seq(seq[:MAX_SEQ_LEN]) ) # "A G S T ..."
    f.close()
    return data

def load_abundance(input_filepath):
    data = []
    f = open( input_filepath )
    for line in f.readlines():
        data.append( line )
    f.close()
    return data


def build_amino_acid_tokenizer(input_filepath, tokenizer_path):
    data = load_amino_acid_data(input_filepath)
    build_tokenizer_from_data(data, tokenizer_path)

def load_mRNA_data(input_filepath):
    data = []
    f = open( input_filepath )
    for line in f.readlines():
        if line[0] == ">": continue
        mRNA = line.strip()[:(MAX_SEQ_LEN*3)]
        data.append( split_word_mRNA(mRNA) ) # "ATG GTC GTT GGA ..."
    f.close()
    return data

def build_mRNA_tokenizer(input_filepath, tokenizer_path):
    data = load_mRNA_data(input_filepath)
    build_tokenizer_from_data(data, tokenizer_path)

def read_mRNA(filepath):
    fa_dic = {}
    f = open( filepath )
    id = ""
    seq = ""
    for line in f.readlines():
        if line[0] == ">":
            if seq != "":
                fa_dic[id] = seq
            id = line[1:].strip()
            seq = ""
        else:
            seq += line.strip()
    f.close()
    if seq != "":
        fa_dic[id] = seq
    return fa_dic

def translate(mRNA):
    seq = ""
    for i in range(0, len(mRNA),3):
        codon = mRNA[i:i+3]
        seq += codon_dic.get(codon, "*")
    return seq

def similarity(seq1, seq2):
    #assert( len(seq1) == len(seq2) )
    # compare until the length of seq1
    match_cnt = 0.0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            match_cnt += 1
    return match_cnt / len(seq1)

def save_mRNA_AA( mRNA_dic, aa_filepath, mRNA_filepath ):
    fo = open( aa_filepath, "w" )
    fo2 = open( mRNA_filepath, "w" )
    for id in mRNA_dic:
        mRNA = mRNA_dic[id]
        seq = translate(mRNA)
        fo.write( ">%s\n%s\n" % (id, seq))
        fo2.write( ">%s\n%s\n" % (id, mRNA))
    fo.close()
    fo2.close()

def load_aminoacid_mRNA_data(aa_input_filepath, mRNA_input_filepath):
    from datasets import Dataset
    aa_data = load_amino_acid_data(aa_input_filepath)
    mRNA_data = load_mRNA_data(mRNA_input_filepath)
    my_list = []
    for i in range(len(aa_data)):
        my_list.append( {"translation":{"aa":aa_data[i],"mRNA":mRNA_data[i]}} )
    dataset = Dataset.from_list(my_list, split='train')
    return dataset

def single_prediction(config, seq, mRNA):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("Using device:", device)

    my_list = []

    words_seq = split_word_amino_acid_seq( seq )
    words_mRNA = split_word_amino_acid_seq( mRNA )
    my_list.append( {"translation":{"aa":words_seq,"mRNA":words_mRNA}} )
    my_list.append( {"translation":{"aa":words_seq,"mRNA":words_mRNA}} )
    from datasets import Dataset
    
    dataset = Dataset.from_list(my_list)

    ## In case tokenizer_src is not ready
    #src_tokenizer_path = Path(config['tokenizer_file'].format("aa"))
    #tokenizer_src = Tokenizer.from_file(str(src_tokenizer_path))
    #trg_tokenizer_path = Path(config['tokenizer_file'].format("mRNA"))
    #tokenizer_tgt = Tokenizer.from_file(str(trg_tokenizer_path))

    model, tokenizer_src, tokenizer_tgt = load_model()

    val_ds = BilingualDataset(dataset, tokenizer_src, tokenizer_tgt, config['lang_src'], config['lang_tgt'], config['seq_len'])
    val_dataloader = DataLoader(val_ds, batch_size=1, shuffle=True)


    if True:
        count = 0
        source_texts = []
        expected = []
        predicted = []
        for batch in val_dataloader:
            print(batch)
            count += 1
            print("ITERATION COUNT =", count)
            encoder_input = batch["encoder_input"].to(device) # (b, seq_len)
            encoder_mask = batch["encoder_mask"].to(device) # (b, 1, 1, seq_len)

            # check that the batch size is 1
            assert encoder_input.size(
                0) == 1, "Batch size must be 1 for validation"

            model_out = greedy_decode(model, encoder_input, encoder_mask, tokenizer_src, tokenizer_tgt, config['seq_len'], device)

            source_text = batch["src_text"][0]
            #target_text = batch["tgt_text"][0]
            model_out_text = tokenizer_tgt.decode(model_out.detach().cpu().numpy())
            model_out_text += "[EOS]"

            source_texts.append(source_text)
            #expected.append(target_text)
            predicted.append(model_out_text)
            
            # Print the source, target and model output
            print(f"{f'SOURCE: ':>12}{source_text}")
            #print(f"{f'TARGET: ':>12}{target_text}")
            print(f"{f'PREDICTED: ':>12}{model_out_text}")

            #mRNA_original = target_text.replace(" ", "")
            mRNA_predicted = model_out_text.replace(" ", "")
            aa_predicted = translate( mRNA_predicted )    
            
            print(f"{f'TARGET AAs: ':>12}{seq}")
            print(f"{f'PREDICTED AAs: ':>12}{aa_predicted}")
            print(f"{f'SIMILARITY SCORE: ':>12}{similarity( seq, aa_predicted )}")
            


def get_ds(config):

    ds_raw = load_aminoacid_mRNA_data("data/chlamydomonas_aa.fa", "data/chlamydomonas_mRNA.fa")

    # Build tokenizers
    tokenizer_src = get_or_build_tokenizer(config, ds_raw, config['lang_src'])

    # test
    #print( config['lang_src'] )
    #i = 0
    #for item in ds_raw:
    #    if ( i > 5 ): break
    #    print( item['translation'][config['lang_src']] )
    #    i += 1
    tokenizer_tgt = get_or_build_tokenizer(config, ds_raw, config['lang_tgt'])

    # test
    #print( config['lang_tgt'] )
    #i = 0
    #for item in ds_raw:
    #    if ( i > 5 ): break
    #    print( item['translation'][config['lang_tgt']] )
    #    i += 1
    
    # Keep 90% for training, 10% for validation
    train_ds_size = int(0.9 * len(ds_raw))
    val_ds_size = len(ds_raw) - train_ds_size
    train_ds_raw, val_ds_raw = random_split(ds_raw, [train_ds_size, val_ds_size])

    train_ds = BilingualDataset(train_ds_raw, tokenizer_src, tokenizer_tgt, config['lang_src'], config['lang_tgt'], config['seq_len'])
    val_ds = BilingualDataset(val_ds_raw, tokenizer_src, tokenizer_tgt, config['lang_src'], config['lang_tgt'], config['seq_len'])

    # Find the maximum length of each sentence in the source and target sentence
    max_len_src = 0
    max_len_tgt = 0

    for item in ds_raw:
        src_ids = tokenizer_src.encode(item['translation'][config['lang_src']]).ids
        tgt_ids = tokenizer_tgt.encode(item['translation'][config['lang_tgt']]).ids
        max_len_src = max(max_len_src, len(src_ids))
        max_len_tgt = max(max_len_tgt, len(tgt_ids))

    print(f'Max length of source sentence: {max_len_src}')
    print(f'Max length of target sentence: {max_len_tgt}')
    
    train_dataloader = DataLoader(train_ds, batch_size=config['batch_size'], shuffle=True)
    val_dataloader = DataLoader(val_ds, batch_size=1, shuffle=True)

    return train_dataloader, val_dataloader, tokenizer_src, tokenizer_tgt

def get_model(config, vocab_src_len, vocab_tgt_len):
    #print (vocab_src_len) # 15698  # 20 AAs
    #print (vocab_tgt_len) # 22463  # 64 Codons
    model = build_transformer(vocab_src_len, vocab_tgt_len, config["seq_len"], config['seq_len'], d_model=config['d_model'], N=config['N'], d_ff=config['d_ff'])
    return model

def load_model():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("Using device:", device)

    ## In case tokenizer_src is not ready
    src_tokenizer_path = Path(config['tokenizer_file'].format("aa"))
    tokenizer_src = Tokenizer.from_file(str(src_tokenizer_path))
    trg_tokenizer_path = Path(config['tokenizer_file'].format("mRNA"))
    tokenizer_tgt = Tokenizer.from_file(str(trg_tokenizer_path))

    model = get_model(config, tokenizer_src.get_vocab_size(), tokenizer_tgt.get_vocab_size()).to(device)

    model_number = 19
    model_filename = get_weights_file_path(config, model_number)
    print(f'Preloading model {model_filename}')
    state = torch.load(model_filename)
    model.load_state_dict(state['model_state_dict'])
    initial_epoch = state['epoch'] + 1
    optimizer = torch.optim.Adam(model.parameters(), lr=config['lr'], eps=1e-9) ## Not used
    optimizer.load_state_dict(state['optimizer_state_dict']) ## Not used
    global_step = state['global_step'] ## Not used

    return model, tokenizer_src, tokenizer_tgt


def train_model(config):
    # Define the device
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("Using device:", device)

    # Make sure the weights folder exists
    Path(config['model_folder']).mkdir(parents=True, exist_ok=True)

    train_dataloader, val_dataloader, tokenizer_src, tokenizer_tgt = get_ds(config)
    model = get_model(config, tokenizer_src.get_vocab_size(), tokenizer_tgt.get_vocab_size()).to(device)
    # Tensorboard
    writer = SummaryWriter(config['experiment_name'])

    optimizer = torch.optim.Adam(model.parameters(), lr=config['lr'], eps=1e-9)

    # If the user specified a model to preload before training, load it
    initial_epoch = 0
    global_step = 0
    if config['preload']:
        model_filename = get_weights_file_path(config, config['preload'])
        print(f'Preloading model {model_filename}')
        state = torch.load(model_filename)
        model.load_state_dict(state['model_state_dict'])
        initial_epoch = state['epoch'] + 1
        optimizer.load_state_dict(state['optimizer_state_dict'])
        global_step = state['global_step']

        # Run validation at the end of every epoch
        batch_iterator = tqdm(train_dataloader, desc=f"Preload Validation")
        run_validation(model, val_dataloader, tokenizer_src, tokenizer_tgt, config['seq_len'], device, lambda msg: batch_iterator.write(msg), global_step, writer)

    loss_fn = nn.CrossEntropyLoss(ignore_index=tokenizer_src.token_to_id('[PAD]'), label_smoothing=0.1).to(device)

    for epoch in range(initial_epoch, config['num_epochs']):
        model.train()
        batch_iterator = tqdm(train_dataloader, desc=f"Processing Epoch {epoch:02d}")
        for batch in batch_iterator:

            encoder_input = batch['encoder_input'].to(device) # (b, seq_len)
            decoder_input = batch['decoder_input'].to(device) # (B, seq_len)
            encoder_mask = batch['encoder_mask'].to(device) # (B, 1, 1, seq_len)
            decoder_mask = batch['decoder_mask'].to(device) # (B, 1, seq_len, seq_len)

            # Run the tensors through the encoder, decoder and the projection layer
            encoder_output = model.encode(encoder_input, encoder_mask) # (B, seq_len, d_model)
            decoder_output = model.decode(encoder_output, encoder_mask, decoder_input, decoder_mask) # (B, seq_len, d_model)
            proj_output = model.project(decoder_output) # (B, seq_len, vocab_size)

            # Compare the output with the label
            label = batch['label'].to(device) # (B, seq_len)

            # Compute the loss using a simple cross entropy
            loss = loss_fn(proj_output.view(-1, tokenizer_tgt.get_vocab_size()), label.view(-1))
            batch_iterator.set_postfix({"loss": f"{loss.item():6.3f}"})

            # Log the loss
            writer.add_scalar('train loss', loss.item(), global_step)
            writer.flush()

            # Backpropagate the loss
            loss.backward()

            # Update the weights
            optimizer.step()
            optimizer.zero_grad(set_to_none=True)

            global_step += 1

        # Run validation at the end of every epoch
        run_validation(model, val_dataloader, tokenizer_src, tokenizer_tgt, config['seq_len'], device, lambda msg: batch_iterator.write(msg), global_step, writer)

        # Save the model at the end of every epoch
        model_filename = get_weights_file_path(config, f"{epoch:02d}")
        torch.save({
            'epoch': epoch,
            'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'global_step': global_step
        }, model_filename)


def generate_mRNA_AA_data():
    mRNA_dic = read_mRNA(filepath = "./data/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.cds.all.fa")
    save_mRNA_AA( mRNA_dic, "data/chlamydomonas_aa.fa", "data/chlamydomonas_mRNA.fa" )
    build_amino_acid_tokenizer("data/chlamydomonas_aa.fa", "tokenizer_aa.json")
    build_mRNA_tokenizer("data/chlamydomonas_mRNA.fa", "tokenizer_mRNA.json")
    
if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    config = get_config()

    if not Path.exists(Path("tokenizer_aa.json")):    # to make tokenizer, etc..
        generate_mRNA_AA_data()
    
    
    if False: #manual switch for train or single_pred
        train_model(config)
    else:
        #random sequence for test
        seq = "MGQQPGKVLGDQRRPSLPALHFIKGAGKRDSSRHGGPHCNVFVEHEALQRPVASDFEPQGLSEAARWNSKENLLAGPSENDPNLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVNSLEKHSWYHGPVSRNAAEYLLSSGINGSFLVRESESSPGQRSISLRYEGRVYHYRINTASDGKLYVSSESRFNTLAELVHHHSTVADGLITTLHYPAPKRNKPTVYGVSPNYDKWEMERTDITMKHKLGGGQYGEVYEGVWKKYSLTVAVKTLKEDTMEVEEFLKEAAVMKEIKHPNLVQLLGVCTREPPFYIITEFMTYGNLLDYLRECNRQEVNAVVLLYMATQISSAMEYLEKKNFIHRNLAARNCLVGENHLVKVADFGLSRLMTGDTYTAHAGAKFPIKWTAPESLAYNKFSIKSDVWAFGVLLWEIATYGMSPYPGIDLSQVYELLEKDYRMERPEGCPEKVYELMRACWQWNPSDRPSFAEIHQAFETMFQESSISDEVEKELGKENLYFQ*"
        mRNA = ""
        single_prediction(config, seq[:40], mRNA[:120])
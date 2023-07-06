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

def read_mRNA(filepath = "./data/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.cds.all.fa"):
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
            target_text = batch["tgt_text"][0]
            model_out_text = tokenizer_tgt.decode(model_out.detach().cpu().numpy())


            source_texts.append(source_text)
            expected.append(target_text)
            predicted.append(model_out_text)
            
            # Print the source, target and model output
            print(f"{f'SOURCE: ':>12}{source_text}")
            print(f"{f'TARGET: ':>12}{target_text}")
            print(f"{f'PREDICTED: ':>12}{model_out_text}")

            mRNA_original = target_text.replace(" ", "")
            mRNA_predicted = model_out_text.replace(" ", "")
            aa_original = translate( mRNA_original )
            aa_predicted = translate( mRNA_predicted )    
            
            print(f"{f'TARGET AAs: ':>12}{aa_original}")
            print(f"{f'PREDICTED AAs: ':>12}{aa_predicted}")
            print(f"{f'SIMILARITY SCORE: ':>12}{similarity( aa_original, aa_predicted )}")
            


def get_ds(config):
    # It only has the train split, so we divide it overselves
    #ds_raw = load_dataset('opus_books', f"{config['lang_src']}-{config['lang_tgt']}", split='train')
    ds_raw_en_it = load_dataset('opus_books', f"en-it", split='train')
    #print(ds_raw_en_it.info)
    '''DatasetInfo(description='This is a collection of copyright free books aligned by Andras Farkas, which are available from http://www.farkastranslations.com/bilingual_books.php\nNote that the texts are rather dated due to copyright issues and that some of them are manually reviewed (check the meta-data at the top of the corpus files in XML). The source is multilingually aligned, which is available from http://www.farkastranslations.com/bilingual_books.php. In OPUS, the alignment is formally bilingual but the multilingual alignment can be recovered from the XCES sentence alignment files. Note also that the alignment units from the original source may include multi-sentence paragraphs, which are split and sentence-aligned in OPUS.\nAll texts are freely available for personal, educational and research use. Commercial use (e.g. reselling as parallel books) and mass redistribution without explicit permission are not granted. Please acknowledge the source when using the data!\n\n16 languages, 64 bitexts\ntotal number of files: 158\ntotal number of tokens: 19.50M\ntotal number of sentence fragments: 0.91M\n', citation="@InProceedings{TIEDEMANN12.463,\n  author = {J�rg Tiedemann},\n  title = {Parallel Data, Tools and Interfaces in OPUS},\n  
    booktitle = {Proceedings of the Eight International Conference on Language Resources and Evaluation (LREC'12)},\n  year = {2012},\n  month = {may},\n  date = {23-25},\n  address = {Istanbul, Turkey},\n  editor = {Nicoletta Calzolari (Conference Chair) and Khalid Choukri and Thierry Declerck and Mehmet Ugur Dogan and Bente Maegaard and Joseph Mariani and Jan Odijk and Stelios Piperidis},\n  publisher = {European Language Resources Association (ELRA)},\n  isbn = {978-2-9517408-7-7},\n  language = {english}\n }\n", 
    homepage='http://opus.nlpl.eu/Books.php', license='', 
    features={'id': Value(dtype='string', id=None), 'translation': Translation(languages=['en', 'it'], id=None)}, 
    post_processed=None, supervised_keys=None, task_templates=None, builder_name='opus_books', config_name='en-it', version=1.0.0, 
    splits={'train': SplitInfo(name='train', num_bytes=8993755, num_examples=32332, shard_lengths=None, dataset_name='opus_books')}, 
    download_checksums={'https://object.pouta.csc.fi/OPUS-Books/v1/moses/en-it.txt.zip': {'num_bytes': 3295251, 'checksum': None}},
    download_size=3295251, post_processing_size=None, dataset_size=8993755, size_in_bytes=12289006)
    '''

    ds_raw = load_aminoacid_mRNA_data("data/chlamydomonas_aa.fa", "data/chlamydomonas_mRNA.fa")
    #print(ds_raw.info)
    '''DatasetInfo(description='', citation='', homepage='', license='', 
    features={'translation': {'aa': Value(dtype='string', id=None), 'mRNA': Value(dtype='string', id=None)}}, 
    post_processed=None, supervised_keys=None, task_templates=None, builder_name=None, config_name=None, version=None, 
    splits=None, 
    download_checksums=None, download_size=None, post_processing_size=None, dataset_size=None, size_in_bytes=None)
    '''

    # Build tokenizers
    tokenizer_src = get_or_build_tokenizer(config, ds_raw, config['lang_src'])

    # test
    #print( config['lang_src'] )
    #i = 0
    #for item in ds_raw:
    #    if ( i > 5 ): break
    #    print( item['translation'][config['lang_src']] )
    #    i += 1
    '''
    Source: Project Gutenberg
    Jane Eyre
    Charlotte Bronte
    CHAPTER I
    There was no possibility of taking a walk that day.
    We had been wandering, indeed, in the leafless shrubbery an hour in the morning; but since dinner (Mrs. Reed, when there was no company, dined early) the cold winter wind had brought with it clouds so sombre, and a rain so penetrating, that further out-door exercise was now out of the question.
    I was glad of it: I never liked long walks, especially on chilly afternoons: dreadful to me was the coming home in the raw twilight, with nipped fingers and toes, and a heart saddened by the chidings of Bessie, the nurse, and humbled by the consciousness of my physical inferiority to Eliza, John, and Georgiana Reed.
    The said Eliza, John, and Georgiana were now clustered round their mama in the drawing-room: she lay reclined on a sofa by the fireside, and with her darlings about her (for the time neither quarrelling nor crying) looked perfectly happy.
    Me, she had dispensed from joining the group; saying, "She regretted to be under the necessity of keeping me at a distance; but that until she heard from Bessie, and could discover by her own observation, that I was endeavouring in good earnest to acquire a more sociable and childlike disposition, a more attractive and sprightly manner--something lighter, franker, more natural, as it were--she really must exclude me from privileges intended only for contented, happy, little children."
    "What does Bessie say I have done?" I asked.
    "Jane, I don't like cavillers or questioners; besides, there is something truly forbidding in a child taking up her elders in that manner.
    Be seated somewhere; and until you can speak pleasantly, remain silent."
    A breakfast-room adjoined the drawing-room, I slipped in there.
    It contained a bookcase: I soon possessed myself of a volume, taking care that it should be one stored with pictures.
    I mounted into the window- seat: gathering up my feet, I sat cross-legged, like a Turk; and, having drawn the red moreen curtain nearly close, I was shrined in double retirement.
    Folds of scarlet drapery shut in my view to the right hand; to the left were the clear panes of glass, protecting, but not separating me from the drear November day.
    At intervals, while turning over the leaves of my book, I studied the aspect of that winter afternoon. Afar, it offered a pale blank of mist and cloud; near a scene of wet lawn and storm-beat shrub, with ceaseless rain sweeping away wildly before a long and lamentable blast.
    I returned to my book--Bewick's History of British Birds: the letterpress thereof I cared little for, generally speaking; and yet there were certain introductory pages that, child as I was, I could not pass quite as a blank.
    They were those which treat of the haunts of sea-fowl; of "the solitary rocks and promontories" by them only inhabited; of the coast of Norway, studded with isles from its southern extremity, the Lindeness, or Naze, to the North Cape-- "Where the Northern Ocean, in vast whirls, Boils round the naked, melancholy isles Of farthest Thule; and the Atlantic surge Pours in among the stormy Hebrides."
    Nor could I pass unnoticed the suggestion of the bleak shores of Lapland, Siberia, Spitzbergen, Nova Zembla, Iceland, Greenland, with "the vast sweep of the Arctic Zone, and those forlorn regions of dreary space,--that reservoir of frost and snow, where firm fields of ice, the accumulation of centuries of winters, glazed in Alpine heights above heights, surround the pole, and concentre the multiplied rigours of extreme cold."
    Of these death-white realms I formed an idea of my own: shadowy, like all the half-comprehended notions that float dim through children's brains, but strangely impressive.
    '''
    tokenizer_tgt = get_or_build_tokenizer(config, ds_raw, config['lang_tgt'])

    # test
    #print( config['lang_tgt'] )
    #i = 0
    #for item in ds_raw:
    #    if ( i > 5 ): break
    #    print( item['translation'][config['lang_tgt']] )
    #    i += 1
    '''
    Source: www.liberliber.it/Audiobook available here
    Jane Eyre
    Charlotte Brontë
    PARTE PRIMA
    I. In quel giorno era impossibile passeggiare.
    La mattina avevamo errato per un'ora nel boschetto spogliato di foglie, ma dopo pranzo (quando non vi erano invitati, la signora Reed desinava presto), il vento gelato d'inverno aveva portato seco nubi così scure e una pioggia così penetrante, che non si poteva pensare a nessuna escursione.
    Ne ero contenta. Non mi sono mai piaciute le lunghe passeggiate, sopra tutto col freddo, ed era cosa penosa per me di tornar di notte con le mani e i piedi gelati, col cuore amareggiato dalle sgridate di Bessie, la bambinaia, e con lo spirito abbattuto dalla coscienza della mia inferiorità fisica di fronte a Eliza, a John e a Georgiana Reed.
    Eliza, John e Georgiana erano aggruppati in salotto attorno alla loro mamma; questa, sdraiata sul sofà accanto al fuoco e circondata dai suoi bambini, che in quel momento non questionavano fra loro né piangevano, pareva perfettamente felice.
    Ella mi aveva proibito di unirmi al loro gruppo, dicendo che deplorava la necessità in cui trovavasi di tenermi così lontana, ma che fino al momento in cui Bessie non guarentirebbe che mi studiavo di acquistare un carattere più socievole e più infantile, maniere più cortesi e qualcosa di più radioso, di più aperto, di più sincero, non poteva concedermi gli stessi privilegi che ai bambini allegri e soddisfatti.
    — Che cosa vi ha detto Bessie di nuovo sul conto mio? — domandai.
    — Jane, non mi piace di essere interrogata. Sta male, del resto, che una bimba tratti così i suoi superiori.
    Sedetevi in qualche posto e state buona fino a quando non saprete parlare ragionevolmente.
    Una piccola sala da pranzo metteva nel salotto, andai in quella pian piano.
    Vi era una biblioteca e io m'impossessai di un libro, cercando che fosse ornato d'incisioni.
    Mi collocai allora nel vano di una finestra, sedendomi sui piedi come i turchi, e tirando la tenda di damasco rosso, mi trovai rinchiusa in un doppio ritiro.
    Le larghe pieghe della cortina scarlatta mi nascondevano tutto ciò che era alla mia destra: alla mia sinistra una invetriata mi proteggeva, ma non mi separava da una triste giornata di novembre.   
    Di tanto in tanto, sfogliando il libro, gettavo un'occhiata al difuori e studiavo l'aspetto di quella serata d'inverno; in lontananza si scorgeva una pallida striscia di nebbia con nuvole, più vicino alberi bagnati, piante sradicate dal temporale e, infine, una pioggia incessante, che lunghe e lamentevoli ventate respingevano sibilando.
    Tornavo allora al mio libro; era La storia degli uccelli dell'Inghilterra, scritta da Berwich. In generale non mi occupavo del testo, nondimeno c'erano delle pagine d'introduzione che non potevo lasciar passare inosservate, malgrado la mia gioventù.
    Esse parlavano di quei rifugi degli uccelli marini, di quei promontori, di quelle rocce deserte abitate da essi soli, di quelle coste della Norvegia sparse d'isole dalla più meridionale punta al capo più nordico, là dove "l'Oceano Polare mugge in vasti turbini attorno all'isola arida e malinconica di Tule, là ove il mare Atlantico si precipita in mezzo alle Ebridi tempestose."
    Non potevo neppure saltare la descrizione di quei pallidi paesaggi della Siberia, dello Spitzberg, della Nuova-Zembla, dell'Islanda, della verde Finlandia! Ero assorta nel pensiero di quella solitudine della zona artica, di quelle immense regioni abbandonate, di quei serbatoi di ghiaccio, ove i campi di neve accumulati durante gli inverni di molti secoli, ammucchiano montagne su montagne per circondare il polo e vi concentrano tutti i rigori del freddo più intenso.
    Mi ero formata un'idea tutta mia di quei regni pallidi come la morte, idea vaga, come sono tutte le cose capite per metà, che fluttuano nella testa dei bimbi; ma quella che mi figuravo produceva in me uno strano effetto.
    '''
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
    mRNA_dic = read_mRNA()
    save_mRNA_AA( mRNA_dic, "data/chlamydomonas_aa.fa", "data/chlamydomonas_mRNA.fa" )
    build_amino_acid_tokenizer("data/chlamydomonas_aa.fa", "tokenizer_aa.json")
    build_mRNA_tokenizer("data/chlamydomonas_mRNA.fa", "tokenizer_mRNA.json")
    
if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    config = get_config()

    if not Path.exists(Path("tokenizer_aa.json")):    # to make tokenizer, etc..
        generate_mRNA_AA_data()
    
    if False:
        train_model(config)
    else:
        seq = "MVNVPKTKRAFCKGCKKHMMMKVTQYKTGKASLYAQGKRRYDRKQSGYGGQTKPVFHKKAKTTKKIVLRMQCQECKQTCMKGLKRCKHFEIGGDKKKGN*"
        mRNA = "ATGGTGAACGTTCCTAAGACCAAGCGGGCGTTCTGCAAGGGGTGCAAGAAGCACATGATGATGAAGGTCACCCAGTACAAGACTGGCAAGGCCTCCCTCTACGCGCAGGGCAAGCGCCGCTACGACCGCAAGCAGTCGGGTTACGGTGGTCAGACCAAGCCCGTCTTCCACAAGAAGGCCAAGACCACCAAGAAGATCGTGCTGCGCATGCAGTGCCAAGAGTGCAAGCAGACCTGCATGAAGGGCCTGAAGCGCTGCAAGCACTTCGAGATCGGTGGTGACAAGAAGAAGGGCAACTAA"
        single_prediction(config, seq[:40], mRNA[:120])
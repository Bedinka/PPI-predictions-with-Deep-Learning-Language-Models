import argparse
import time
import torch
from Models import get_model
from Process import *
import torch.nn.functional as F
from Optim import CosineWithRestarts
from Batch import create_masks
import pdb
import dill as pickle
import argparse
from Models import get_model
from Beam import beam_search2
from nltk.corpus import wordnet
from torch.autograd import Variable
import re
import os


def multiple_replace(dict, text):
  # Create a regular expression  from the dictionary keys
  regex = re.compile("(%s)" % "|".join(map(re.escape, dict.keys())))

  # For each match, look-up corresponding value in dictionary
  return regex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], text) 

def translate_sentence(sentence, model, opt, SRC, TRG):
    model.eval()
    indexed = []

    for tok in sentence:
        if SRC.vocab.stoi[tok] != 0 or opt.floyd == True:
            indexed.append(SRC.vocab.stoi[tok])

    sentence = Variable(torch.LongTensor([indexed]))
    if opt.device == 0:
        sentence = sentence.cuda()

    sentence = beam_search2(sentence, model, SRC, TRG, opt)

    return multiple_replace({' ?' : '?',' !':'!',' .':'.','\' ':'\'',' ,':','}, sentence)

def translate(opt, model, SRC, TRG):
    #sentences = opt.text.lower().split('.')
    sentences = opt.text.split()
    translated = []

    for sentence in sentences:
        #translated.append(translate_sentence(sentence + '.', model, opt, SRC, TRG).capitalize())
        translated.append(translate_sentence(sentence, model, opt, SRC, TRG).capitalize())

    return (' '.join(translated))


def translate_file(file_path, model, SRC, TRG, opt):
    try:
        with open(file_path, encoding='utf-8') as file:
            content = file.read()
    except:
        print("Error opening or reading the file.")
        return

    lines = content.split('\n')  # Split content into separate lines

    translated_lines = []  # Store the translated lines

    total_lines = len(lines)
    processed_lines = 0

    for line in lines:
        if line:  # Skip empty lines
            translated_line = translate_sentence(line, model, opt, SRC, TRG).capitalize()
            formatted_translation = translated_line.split('<eos>')[0].replace(' ', '')[5:] + '\n'
            translated_lines.append(formatted_translation)

            processed_lines += 1
            print(f"Processed lines: {processed_lines}/{total_lines}")

    # Join the translated lines into a single string
    translated_content = ''.join(translated_lines)

    # Create a new file path for saving the translation
    file_name, file_extension = os.path.splitext(file_path)
    output_file_path = file_name + "_translated" + file_extension

    try:
        with open(output_file_path, 'w', encoding='utf-8') as output_file:
            output_file.write(translated_content)
    except:
        print("Error writing to the output file.")
        return

    print("Translation saved to:", output_file_path)


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-load_weights', required=True)
    parser.add_argument('-k', type=int, default=3)
    parser.add_argument('-max_len', type=int, default=1500)
    parser.add_argument('-d_model', type=int, default=512)
    parser.add_argument('-n_layers', type=int, default=6)
    parser.add_argument('-heads', type=int, default=8)
    parser.add_argument('-dropout', type=int, default=0.1)
    parser.add_argument('-no_cuda', action='store_true')
    parser.add_argument('-floyd', action='store_true')
    
    opt = parser.parse_args()

    opt.device = 0 if opt.no_cuda is False else -1
 
    assert opt.k > 0
    assert opt.max_len > 10

    SRC, TRG = create_fields(opt)
    #print(len(SRC.vocab), len(TRG.vocab))
    #print(SRC.vocab)
    #print(dir(SRC.vocab))
    #print(SRC.vocab.__dict__)
    #print(SRC.vocab.stoi["H"])
    
    model = get_model(opt, len(SRC.vocab), len(TRG.vocab))
    
    while True:
        opt.text =input("Enter a sentence to translate (type 'f' to load from file, or 'q' to quit):\n")
        if opt.text=="q":
            break
        if opt.text=='f':
            file_path = input("Enter the file path to translate:\n")
            translate_file(file_path, model, SRC, TRG, opt)
            continue

        print("opt.text :", opt.text)
        phrase = translate(opt, model, SRC, TRG)
        print('> '+ phrase.split('<eos>')[0].replace(' ', '')[5:] + '\n')

if __name__ == '__main__':
    main()

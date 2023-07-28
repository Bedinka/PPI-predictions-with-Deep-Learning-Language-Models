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
from Models import get_model
from Beam import beam_search2
from nltk.corpus import wordnet
from torch.autograd import Variable
import re
import os
from Bio import AlignIO
from Bio.Align import PairwiseAligner
import statistics
import numpy as np

def multiple_replace(dict, text):
    # Create a regular expression from the dictionary keys
    regex = re.compile("(%s)" % "|".join(map(re.escape, dict.keys())))

    # For each match, look up the corresponding value in the dictionary
    return regex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], text)


def calculate_alignment_scores(predictions, originals):
    alignment_scores = []
    noncoincidences = []

    for prediction, original in zip(predictions, originals):
        score, noncoincidence = calculate_match_score(prediction, original)
        alignment_scores.append(score)
        noncoincidences.append(noncoincidence)

    return alignment_scores, noncoincidences


def calculate_match_score(prediction, original):
    match_count = sum(1 for a, b in zip(prediction, original) if a == b)
    score = match_count / len(original) * 100.0
    noncoincidence = len(original) - match_count
    return score, noncoincidence


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
    translated_sentence = multiple_replace({' ?': '?', ' !': '!', ' .': '.', '\' ': '\'', ' ,': ','}, sentence).upper()

    return translated_sentence


def translate(opt, model, SRC, TRG, original_sequences):
    sentences = opt.text.split()
    translated = []

    for sentence, original in zip(sentences, original_sequences):
        translated.append(translate_sentence(sentence, model, opt, SRC, TRG).capitalize())

    return ' '.join(translated)


def load_original_sequences(file_path):
    try:
        with open(file_path, encoding='utf-8') as file:
            content = file.read()
    except:
        print("Error opening or reading original file.")
        return []

    lines = content.split('\n')  # Split content into separate lines
    sequences = []
    current_sequence = ""

    for line in lines:
        if line and not line.startswith('>'):  # Process non-empty lines that are not ID lines
            current_sequence += line.strip()

        elif line:  # Process lines starting with ">"
            if current_sequence:
                sequences.append(current_sequence)
            current_sequence = ""

    if current_sequence:
        sequences.append(current_sequence)

    return sequences


def translate_file(file_path, original_file_path, model, SRC, TRG, opt):
    original_sequences = load_original_sequences(original_file_path)

    try:
        with open(file_path, encoding='utf-8') as file:
            content = file.read()
    except:
        print("Error opening or reading the file.")
        return

    lines = content.split('\n')  # Split content into separate lines

    translated_lines = []  # Store the translated lines
    ids = []  # Store the corresponding IDs

    total_lines = len(lines)
    processed_lines = 0

    current_id = None
    current_translation = ""
    current_original = ""
    predictions = []

    for line in lines:
        if line and line.startswith('>'):  # Process lines starting with ">"
            if current_id is not None:  # Store the previous ID and translation
                ids.append(current_id)
                translated_lines.append(current_translation.strip())
                predictions.append(current_translation.strip())

            current_id = line.strip()[1:]  # Store the new ID
            current_translation = ""  # Reset the translation
            current_original = ""  # Reset the original sequence

        elif line:  # Process non-empty lines that are not ID lines
            translated_line = translate_sentence(line, model, opt, SRC, TRG).capitalize()
            formatted_translation = translated_line.split('<eos>')[0].replace(' ', '')[5:].upper()
            current_translation += formatted_translation
            current_original += line.strip()

        processed_lines += 1
        print(f"Processed lines: {processed_lines}/{total_lines}")

    # Store the last ID and translation
    if current_id is not None:
        ids.append(current_id)
        translated_lines.append(current_translation.strip())
        predictions.append(current_translation.strip())

    # Calculate alignment scores, identities, and non-coincidences
    alignment_scores, noncoincidences = calculate_alignment_scores(predictions, original_sequences)
    

    # Join the translated lines into a single string
    translated_content = '\n'.join(translated_lines)

    # Create a new file path for saving the translation
    file_name, file_extension = os.path.splitext(file_path)
    output_file_path = file_name + "_translated" + file_extension

    try:
        with open(output_file_path, 'w', encoding='utf-8') as output_file:
            for id, translation in zip(ids, translated_lines):
                output_file.write(f">{id}\n")
                output_file.write(translation + '\n')
    except:
        print("Error writing to the output file.")
        return

    # Calculate statistics
    alignment_scores_abs = [abs(score) for score in alignment_scores]

    noncoincidences_mean= statistics.mean(noncoincidences) if noncoincidences else 0.0
    alignment_scores_mean = statistics.mean(alignment_scores_abs) if alignment_scores_abs else 0.0
    alignment_scores_std = statistics.stdev(alignment_scores_abs) if alignment_scores_abs else 0.0
    alignment_scores_var = statistics.variance(alignment_scores_abs) if alignment_scores_abs else 0.0

    # Create a new file path for saving the statistics
    stats_file_path = file_name + "_statistics.txt"

    try:
        with open(stats_file_path, 'w', encoding='utf-8') as stats_file:
            stats_file.write(f"Non-coincidences Mean: {noncoincidences_mean}\n")
            stats_file.write(f"Alignment Scores Mean: {alignment_scores_mean}\n")
            stats_file.write(f"Alignment Scores Standard Deviation: {alignment_scores_std}\n")
            stats_file.write(f"Alignment Scores Variance: {alignment_scores_var}\n")

            
    except:
        print("Error writing the statistics to the file.")
        return

    print("Translation saved to:", output_file_path)
    print("Statistics saved to:", stats_file_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-load_weights', required=True)
    parser.add_argument('-k', type=int, default=3)
    parser.add_argument('-max_len', type=int, default=1500)
    parser.add_argument('-d_model', type=int, default=512)
    parser.add_argument('-n_layers', type=int, default=6)
    parser.add_argument('-heads', type=int, default=8)
    parser.add_argument('-dropout', type=float, default=0.1)
    parser.add_argument('-no_cuda', action='store_true')
    parser.add_argument('-floyd', action='store_true')
    parser.add_argument('-original_file', required=True)

    opt = parser.parse_args()

    opt.device = 0 if opt.no_cuda is False else -1

    assert opt.k > 0
    assert opt.max_len > 10

    SRC, TRG = create_fields(opt)
    model = get_model(opt, len(SRC.vocab), len(TRG.vocab))

    while True:
        opt.text = input("Enter a sentence to translate (type 'f' to load from file, or 'q' to quit):\n")
        if opt.text == "q":
            break
        if opt.text == 'f':
            file_path = input("Enter the file path to translate:\n")
            translate_file(file_path, opt.original_file, model, SRC, TRG, opt)
            continue

        print("opt.text :", opt.text)
        phrase = translate(opt, model, SRC, TRG, [opt.text])
        print('> ' + phrase.split('<eos>')[0].replace(' ', '')[5:].upper() + '\n')


if __name__ == '__main__':
    main()

import argparse
import os
import textwrap
import codon_optimization_score as cs
from train import split_word_amino_acid_seq, similarity, load_model, greedy_decode, generate_mRNA_AA_data, single_prediction, threeway_split
from config import get_config
import warnings
from pathlib import Path
import sys
import torch
from dataset import BilingualDataset
from torch.utils.data import DataLoader
import textwrap
import shap
import pandas as pd
import numpy as np
import shap
""" 
def get_prediction(input_sequence):
    # Convert input_sequence to the format expected by your model
    
    # Make the prediction using your transformer model
    prediction = model.predict(input_sequence)

    # Compute SHAP values
    explainer = shap.Explainer(model)  # Create an explainer object for your model
    shap_values = explainer.shap_values(input_sequence)

    # Visualize the SHAP values or use them for further analysis
    shap.summary_plot(shap_values, input_sequence)  # Example: create a summary plot

    return prediction
 """
def predict_sequence(config, seq, model_number, weights="", basename=""):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("Using device:", device)

    if len(seq) > 100:
        pred = []
        split_list = textwrap.wrap(seq, 100)
        for i in split_list:
            mRNA_predicted, prob_list = single_prediction(config, i, model_number, weights=weights, basename=basename)
            pred.append(mRNA_predicted)
        mRNA_predicted = "".join(pred)
    else:
        mRNA_predicted, prob_list = single_prediction(config, seq, model_number, weights=weights, basename=basename)

    """ model = load_model(model_number)[0]
    # Compute SHAP values
    explainer = shap.Explainer(model)  # Create an explainer object for your model
    shap_values = explainer.shap_values(seq)
    shap.summary_plot(shap_values, seq) """  # Example: create a summary plot
    
    return mRNA_predicted, prob_list

def single_prediction(config, seq, model_number, weights = "", basename = ""):



    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("Using device:", device)

    

    my_list = []

    words_seq = split_word_amino_acid_seq( seq )
    words_mRNA = "[SOS]"
    my_list.append( {"translation":{"aa":words_seq,"mRNA":words_mRNA}} )
    # my_list.append( {"translation":{"aa":words_seq,"mRNA":words_mRNA}} )
    from datasets import Dataset
    
    dataset = Dataset.from_list(my_list)

    ## In case tokenizer_src is not ready
    #src_tokenizer_path = Path(config['tokenizer_file'].format("aa"))
    #tokenizer_src = Tokenizer.from_file(str(src_tokenizer_path))
    #trg_tokenizer_path = Path(config['tokenizer_file'].format("mRNA"))
    #tokenizer_tgt = Tokenizer.from_file(str(trg_tokenizer_path))
    if weights and basename and model_number:
        config["model_folder"] = weights
        config["model_basename"] = basename
        model, tokenizer_src, tokenizer_tgt = load_model(model_number, weights, basename)
    else:
        model, tokenizer_src, tokenizer_tgt = load_model(model_number)
    model.eval()

    # explainer = shap.Explainer(model, tokenizer_src)

    val_ds = BilingualDataset(dataset, tokenizer_src, tokenizer_tgt, config['lang_src'], config['lang_tgt'], config['seq_len'])
    # print(val_ds.ds)

    val_dataloader = DataLoader(val_ds, batch_size=1, shuffle=True)
    # print(val_dataloader.dataset)

    # shap_values = explainer.shap_values(pd.DataFrame(val_dataloader.dataset.ds.src_texts))
    # shap.plots.text(shap_values)

    count = 0
    source_texts = []
    predicted = []
    for batch in val_dataloader:
        # print(batch)
        count += 1
        #print("ITERATION COUNT =", count)
        encoder_input = batch["encoder_input"].to(device) # (b, seq_len)
        encoder_mask = batch["encoder_mask"].to(device) # (b, 1, 1, seq_len)

        # check that the batch size is 1
        assert encoder_input.size(
            0) == 1, "Batch size must be 1 for validation"
        model_out, prob_list = greedy_decode(model, encoder_input, encoder_mask, tokenizer_src, tokenizer_tgt, config['seq_len'], device)
        prob_list = [i.tolist() for i in prob_list]

        source_text = batch["src_text"][0]
        
        model_out_text = tokenizer_tgt.decode(model_out.detach().cpu().numpy()).replace(" ", "")
        model_out_text = model_out_text[:len(seq)*3]
        # model_out_text += "[EOS]"

        source_texts.append(source_text)
        predicted.append(model_out_text)

        mRNA_predicted = model_out_text.replace(" ", "")

        return mRNA_predicted, prob_list
        
def gc_content(mRNA):
    GC = 0.0
    for nt in mRNA:
        if nt in "GC":
            GC += 1
    GC_Content = GC / float( len(mRNA) )
    return GC_Content

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    # parser = argparse.ArgumentParser()
    # parser.add_argument('-s', '--split', action='store_true', help='perform the split behavior')
    # args = parser.parse_args()
    config = get_config()

    if not Path.exists(Path("tokenizer_aa.json")):    # to make tokenizer, etc..
        generate_mRNA_AA_data()

    if len(sys.argv) < 2:
        print("Usage: \n python predict.py <amino acid sequence> \n or \n python predict.py <file name> \n \n Loading example sequence...")
        seq = "MPEGFIKAGQRPSLSGTPLVSANQGVTGMPVSAFTVILSKAYPAIGTPIPFDKILYNRQQHYDPRTGIFTCQIPGIYYFSYHVHVKGTHVWVGLYKNGTPVMYTYDEYTKGYLDQASGSAIIDLTENDQVWLQLPNAESNGLYSSEYVHSSFSGFLVAPM"
    else:
        if os.path.isfile(sys.argv[1]):
            seq = ""
            with open(sys.argv[1], "r") as f:
                for line in f.readlines():
                    if line.startswith(">"):
                        continue
                    else:
                        seq += line.strip()
        else:
            seq = sys.argv[1]

    # seq variable contains the aa sequence
    
    if seq.endswith("*") == False:
        seq += "*" 

    if len(sys.argv) > 2 and sys.argv[2] == "s":
        print("Splitting sequence in 3 parts and predicting each part separately")
        seq3ws = threeway_split(seq, "aa", 26) # threeway split must output the list with the three partitions


        out_start = predict_sequence(config, seq3ws[0], "04", weights = "weights_start", basename = "tmodel_start")
        out_mid = predict_sequence(config, seq3ws[1], "09", weights = "weights_middle", basename = "tmodel_middle")
        out_end = predict_sequence(config, seq3ws[2], "04", weights = "weights_end", basename = "tmodel_end")

        out = "".join([out_start, out_mid, out_end])
        out2 = predict_sequence(config, seq, 19, weights="weights", basename="tmodel_")
        print('NOsplit mRNA: ' + out2)

    elif len(sys.argv) > 2 and sys.argv[2] == "r": #reverse behavior for check
        print("invert start and termination for backward check")
        seq3ws = threeway_split(seq, "aa", 26) # threeway split output the list with the three partitions


        out_start = predict_sequence(config, seq3ws[0], "04", weights = "weights_end", basename = "tmodel_end")
        out_mid = predict_sequence(config, seq3ws[1], "09", weights = "weights_middle", basename = "tmodel_middle")
        out_end = predict_sequence(config, seq3ws[2], "04", weights = "weights_start", basename = "tmodel_start")

        out = "".join([out_start, out_mid, out_end])
        out2 = predict_sequence(config, seq, 19, weights="weights", basename="tmodel_")
        print('NOsplit mRNA: ' + out2)
            
    

    else:
        out, p_list = predict_sequence(config, seq, 19, weights = "weights", basename = "tmodel")
    
    aa_out = cs.translate(out)
    opt_mRNA = cs.optimize_max(aa_out)

    print("pred mRNA: " + out + '\n')
    print("opti mRNA: " + opt_mRNA + '\n')
    print("pred AA: " + aa_out + '\n')
    print("orig AA: " + seq + '\n')

    print("codon optim score: " + str(cs.codon_optimization_score(out)))
    print("GC content: " + str(gc_content(out)))
    print("mRNA similarity score : " + str(similarity(out, opt_mRNA)))
    print("AA similarity score : " + str(similarity(seq, aa_out)))
    print(p_list)


    mRNA_eval, p_list_eval = predict_sequence(config, aa_out, 19, weights="weights", basename="tmodel_")

    error = []
    for i in p_list:
        for j in i:
            error.append(i[j] - p_list_eval[i])

    error = sum(error)/len(error)
    print(error)

    # convert tensor to list
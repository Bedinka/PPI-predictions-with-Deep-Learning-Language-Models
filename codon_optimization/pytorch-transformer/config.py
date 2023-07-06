from pathlib import Path

def get_config():
    return {
        "batch_size": 4, #8,
        "num_epochs": 20,
        "lr": 10**-4,
        "seq_len": 150,
        "d_model": 32, #64, #128, #256, #512,
        "d_ff": 512, #1024, #2048,
        "N": 4, #6, 
        "lang_src": "aa", #"en",
        "lang_tgt": "mRNA", #"it",
        "model_folder": "weights",
        "model_basename": "tmodel_",
        "preload": 19, #None,
        "tokenizer_file": "tokenizer_{0}.json",
        "experiment_name": "runs/tmodel"
    }

def get_weights_file_path(config, epoch: str):
    model_folder = config["model_folder"]
    model_basename = config["model_basename"]
    model_filename = f"{model_basename}{epoch}.pt"
    return str(Path('.') / model_folder / model_filename)



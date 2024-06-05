from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import os 
import warnings
import glob
import pandas as pd

def pdb_to_fasta_double(pdb_file):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_file)
    sequence_A = ''
    sequence_B = ''
    for model in structure:
        for chain in model:
            if chain.get_id() == 'A':
                for residue in chain:
                    if residue.get_id()[0] == ' ':
                        sequence_A += seq1(residue.get_resname())
            elif chain.get_id() == 'B':
                for residue in chain:
                    if residue.get_id()[0] == ' ':
                        sequence_B += seq1(residue.get_resname())
    return sequence_A, sequence_B

def pdb_to_fasta_single(pdb_file):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parser = PDBParser(QUIET=True)  # QUIET=True to suppress warnings
        structure = parser.get_structure('structure', pdb_file)
        sequence = ''
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[0] == ' ':
                        sequence += seq1(residue.get_resname())
    return sequence

def main(input_tsv, pdb_file_dir, interactom_fasta, single ):
    work_dir = "/home/dina/Documents/PPI_WDLLM/workdir"
    if single == True: 
        df = pd.read_csv(input_tsv, sep='\t')
        with open(interactom_fasta, 'a') as f:
            for index, row in df.iterrows():
                protein_ids = row['Protein ID'].split('_')
                if protein_ids[0]:
                    search_pattern = os.path.join(work_dir, pdb_file_dir, f"*{protein_ids[0]}*.pdb")
                    pdb1_files = glob.glob(search_pattern)
                    if pdb1_files:
                        pdb1 = pdb1_files[0]  # Assuming you want the first match
                        sequence = pdb_to_fasta_single(pdb1)
                        f.write(f'>pdb:{protein_ids[0]}_Chain_A\n')
                        f.write(sequence + '\n')
                if protein_ids[1]:
                    search_pattern = os.path.join(work_dir, pdb_file_dir, f"*{protein_ids[1]}*.pdb")
                    pdb2_files = glob.glob(search_pattern)
                    if pdb2_files:
                        pdb2 = pdb2_files[0]  # Assuming you want the first match
                        sequence = pdb_to_fasta_single(pdb2)
                        f.write(f'>pdb:{protein_ids[1]}_Chain_A\n')
                        f.write(sequence + '\n')
    else: 
        processed_pdb_files = [os.path.join(f"{work_dir}/{pdb_file_dir}", file) for file in os.listdir(f'{work_dir}/{pdb_file_dir}') if file.endswith('.pdb')]
        # Open the output file in append mode ('a')
        with open(interactom_fasta, 'a') as f:
            for pdb_file in processed_pdb_files:
                sequence_A, sequence_B = pdb_to_fasta_double(pdb_file)
                f.write('>pdb:A\n')
                f.write(f'>{os.path.basename(pdb_file)}_Chain_A\n')
                f.write(sequence_A + '\n')
                f.write('>pdb:B\n')
                f.write(f'>{os.path.basename(pdb_file)}_Chain_B\n')
                f.write(sequence_B + '\n')
            
if __name__ == "__main__":
    input_tsv="/home/dina/Documents/PPI_WDLLM/2024-05-30_22-44-22/filtered_file_train__2024-05-30_22-44-22_pipelin_testing.tsv"
    interactom_fasta="/home/dina/Documents/PPI_WDLLM/2024-05-30_22-44-22/interactom_train_2024-05-30_22-44-22_pipelin_testing_ca.fasta"
    pdb_file_dir='TRAIN'
    single = False
    main(input_tsv, pdb_file_dir, interactom_fasta, single)
    
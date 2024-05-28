from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import os 
import warnings

def pdb_to_fasta(pdb_file):
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

def main(pdb,interactom_fasta):
    work_dir = "/home/dina/Documents/PPI_WDLLM/workdir"
    processed_pdb_files = [os.path.join(f"{work_dir}/{pdb}", file) for file in os.listdir(f'{work_dir}/{pdb}') if file.endswith('.pdb')]

    # Open the output file in append mode ('a')
    with open(interactom_fasta, 'a') as f:
        for pdb_file in processed_pdb_files:
            sequence_A, sequence_B = pdb_to_fasta(pdb_file)
            f.write('>pdb:A\n')
            f.write(f'>{os.path.basename(pdb_file)}_Chain_A\n')
            f.write(sequence_A + '\n')
            f.write('>pdb:B\n')
            f.write(f'>{os.path.basename(pdb_file)}_Chain_B\n')
            f.write(sequence_B + '\n')
            
if __name__ == "__main__":
    pdb=''
    interactom_fasta=''
    main(pdb, interactom_fasta)
    
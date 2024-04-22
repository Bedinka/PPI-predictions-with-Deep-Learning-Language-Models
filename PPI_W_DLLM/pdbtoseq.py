from Bio.PDB import PDBParser
import os

# Directory containing PDB files
pdb_directory = "/home/dina/Documents/PPI_WDLLM/workdir/interactions_001"

# Function to extract protein sequences from PDB files
def extract_sequences_from_pdb(pdb_directory):
    sequences = []
    parser = PDBParser()
    for filename in os.listdir(pdb_directory):
        if filename.endswith(".pdb"):
            pdb_path = os.path.join(pdb_directory, filename)
            structure = parser.get_structure(filename, pdb_path)
            for model in structure:
                for chain in model:
                    sequence = ""
                    for residue in chain:
                        if residue.get_id()[0] == " ":
                            sequence += residue.get_resname()
                    sequences.append(sequence)
    return sequences

# Write sequences to a tab-separated file
def write_sequences_to_file(sequences, output_file):
    with open(output_file, "w") as f:
        for seq in sequences:
            f.write(seq + "\n")

# Main function
def main():
    sequences = extract_sequences_from_pdb(pdb_directory)
    output_file = "protein_sequences.tsv"
    write_sequences_to_file(sequences, output_file)
    print(f"All protein sequences have been written to {output_file}")

if __name__ == "__main__":
    main()

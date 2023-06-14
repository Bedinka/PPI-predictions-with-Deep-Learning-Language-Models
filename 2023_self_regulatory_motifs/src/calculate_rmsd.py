import sys
from Bio.PDB import CEAligner, PDBParser


def calculate_rmsd( structure1, structure2):
    
    # Calculate RMSD between two protein structures.
    s1 = structure1
    s2 = structure2
    aligner = CEAligner()
    aligner.set_reference(s1)
    aligner.align(s2)
    return aligner.rms

if __name__ == "__main__":
    PDB1 = sys.argv[1]
    PDB2 = sys.argv[2]
    
    parser = PDBParser()
    structure1 = parser.get_structure("PDB1", PDB1)
    structure2 = parser.get_structure("PDB2", PDB2)
        
    print( calculate_rmsd( structure1, structure2 ) )
    #python3 calculate_rmsd.py AF-P56166-F1-model_v4_94-140_PF01620.19.pdb AF-Q40960-F1-model_v4_7-39_PF01620.19.pdb
    
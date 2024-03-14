from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import tarfile
import os
import numpy as np
import pandas as pd
from scipy import stats
import freesasa
import torch 
import pydssp
import tqdm
import sys


# Extracting protein IDs to a txt file 
def extract_pdb_names_from_tgz(tgz_file, output_file):
    pdb_files = []
    with tarfile.open(tgz_file, "r:gz") as tar:
        pdb_files = [member for member in tar.getmembers() if os.path.splitext(member.name)[1] == '.pdb']
        with open(output_file, 'w') as f:
            for pdb_file in pdb_files:
                pdb_name = os.path.splitext(os.path.basename(pdb_file.name))[0]
                first_part, second_part = pdb_name.split('-')[:2]
                f.write(f"{first_part}\t{second_part}\n")

# Unpacking .tgz file 
def process_tgz_files_in_directory(work_dir):
    processed_pdb_files = []
    for tgz_file in os.listdir(work_dir):
        if tgz_file.endswith(".tgz"):
            output_file = os.path.splitext(tgz_file)[0] + ".txt"
            extract_pdb_names_from_tgz(os.path.join(work_dir, tgz_file), os.path.join(work_dir, output_file))
            print(f"Processed {tgz_file} and saved output to {output_file}")

            tgz_dir = os.path.join(work_dir, os.path.splitext(tgz_file)[0])
            os.makedirs(tgz_dir, exist_ok=True)

            with tarfile.open(os.path.join(work_dir, tgz_file), "r:gz") as tar:
                tar.extractall(path=tgz_dir)
                pdb_files = [os.path.join(tgz_dir, member.name) for member in tar.getmembers() if os.path.splitext(member.name)[1] == '.pdb']
                processed_pdb_files.extend(pdb_files)
            
    return processed_pdb_files

# Splitting PDB files by chains 
def splitPDBbyChain(filename, output_dir):
    chains = {}
    splitted_files = []    
    with open(filename, 'r') as f:
        first_three_lines = [next(f) for _ in range(3)]
        for line in f:
            if line[:4] != "ATOM":
                continue
            chainID = line[21]
            if chainID not in chains:
                chains[chainID] = []
            chains[chainID].append(line)

    basename_ext = os.path.basename(filename)
    basename = os.path.splitext(basename_ext)[0]  
    output_subdir = os.path.join(output_dir, f"{basename}_splitted")
    os.makedirs(output_subdir, exist_ok=True)

    for chainID in chains:
        chain_filename = os.path.join(output_subdir, f"{basename}-{chainID}.pdb")
        with open(chain_filename, "w") as fout:
            fout.writelines(first_three_lines)
            for line in chains[chainID]:
                fout.write(line)
        splitted_files.append(chain_filename)

    return splitted_files, 

class Chain:
    def __init__(self, chainID):
        self.chainID = chainID
        self.residues = {}
        self.residue_indexes = []    
        self.prot_id = []
          

    def addResidue(self, aa, resnum):
        aResidue = Residue( aa, resnum )
        self.residues[resnum] = aResidue
        self.residue_indexes.append(resnum)
        #print(aResidue)
        pass
    
    def addCA(self, aa, resnum, x, y, z):
        aResidue = Residue( aa, resnum, x, y, z)
        self.residues[resnum] = aResidue
        self.residue_indexes.append(resnum)
        #print(aResidue)
        pass
    
    def get_all_atoms(self):
        all_atoms = []
        all_coords = []
        for resnum in self.residue_indexes:
            residue = self.residues[resnum]
            for atom_name, coordinates in residue.atoms.items():
                all_atoms.append(atom_name)
                all_coords.append(coordinates)
        return all_atoms, all_coords
    
class Residue:
    def __init__(self, aa, resnum  ):
        self.aa = aa
        self.resnum = resnum
        self.atoms = {}
        self.rsa_value = None
        pass

    def addAtom(self, atom_name, x, y, z):
        self.atoms[atom_name] = (x,y,z)

    def add_coordinate(self, coordinate):
        self.coordinates.append(coordinate)

    def get_mean_coordinate(self):
        xs = []
        ys = []
        zs = []
        for x, y, z in self.atoms.values():
            xs.append(x)
            ys.append(y)
            zs.append(z)
        return ( np.mean(xs), np.mean(ys), np.mean(zs) )
    
    def __str__(self):
        return self.aa + " - " + str(self.resnum) 
    
    def addRSA(self, sasa_all):
        self.rsa_value = sasa_all
        print(f"Added one rsa to {self.resnum}")
        print(self.rsa_value)
        
    
    
class Matrix:
    def __init__(self, row, col):
        self.matrix = []
        self.row =row
        self.col = col
        self.i=0 
        self.j=0
        self.residues = []

    def submatrixes( dist_mat , size, overlap ): 
        submatrices = create_fixedsize_submatrix(dist_mat, size, overlap)
        return submatrices

    def sub_residuses(self, chains_CA, chainID1, chainID2, submatrices):

        self.sub_res_a = []
        self.sub_res_b = []
        Chain_A = chains_CA[chainID1]
        Chain_B = chains_CA[chainID2]

        """ getting parsepdb chains keys as dictionary keys 
        self.sub_res_a = {}
        self.sub_res_b = {}
        for chain_id in chains_CA.keys():
        if chain_id == chainID1:
            self.sub_res_a[chain_id] = []
        elif chain_id == chainID2:
            self.sub_res_b[chain_id] = []
        """

        for idx, submatrix in enumerate(submatrices):
            print(f"Submatrix {idx+1}:")
            print(submatrix)

            # matrix, i, j
            print(submatrix[0]) #matrix
            print(submatrix[1]) #i
            print(submatrix[2]) #j
            
            i = submatrix[1]
            j = submatrix[2]
            
            start_a = Chain_A.residue_indexes[i]
            start_b = Chain_B.residue_indexes[j]
            sub_res_a = Chain_A.residues[start_a]
            sub_res_b = Chain_A.residues[start_b]
    
        return  sub_res_a, sub_res_b

    def __getitem__(self, key):
        return self.matrix[key]
        
# Parsing the PDB files 
def parsePDB(pdb_file):
    chains = {}
    parts = pdb_file.split('/')
    filename = parts[-1]
    filename_parts = filename.split('-')
    first_part = filename_parts[0]
    second_part = filename_parts[1]
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") == False:
                continue
            #if line[12:16].strip() != "CA":
            #    continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            aa = line[17:20] 
            atom_name = line[12:16].strip()
            resnum = int(line[22:26]) 
            chainID = line[21]
            if chainID not in chains:
                chains[chainID] = Chain(chainID)
                if chainID == 'A':
                    chains[chainID].prot_id = first_part
                    print(chains[chainID].prot_id)
                else :
                    chains[chainID].prot_id = second_part
            if resnum not in chains[chainID].residues:
                chains[chainID].addResidue( aa, resnum)
            chains[chainID].residues[resnum].addAtom( atom_name, x, y, z ) # .addCA(aa, resnum, x, y, z)    
    return chains

# Calculating distance 
def calculate_distance(p1, p2):
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    dz = p1[2] - p2[2]
    distsq = pow(dx, 2) + pow(dy, 2) + pow(dz, 2)
    distance = np.sqrt(distsq)
    return distance

# Finding iteracting fragments by given distance
def findInteractingFragments(chains):
    [_chainA, _chainB] = chains.keys()
    chainA = chains[_chainA]
    chainB = chains[_chainB]
    interactions=[]
    for resnumA in chainA.residue_indexes:
        residueA = chainA.residues[resnumA]
        for resnumB in chainB.residue_indexes:
            residueB = chainB.residues[resnumB]
            distance = calculate_distance(residueA.CA, residueB.CA)
            if distance < 11.0:
                print("Found", residueA, residueB, distance)
                interactions.append((residueA, residueB, distance))
    return interactions

# Finding iteracting residues by given distance
def findInteractingResidues(chains, chainID1, chainID2, distance_matrix_CA):
    interactions = []
    chainA = chains[chainID1]
    chainB = chains[chainID2]
    interaction_matrix = np.zeros((len(chainA.residue_indexes),len(chainB.residue_indexes)))
    for i in range(distance_matrix_CA.shape[0]): # row
        for j in range(distance_matrix_CA.shape[1]): # column
            distance = distance_matrix_CA[i,j]
            if distance < 11.0:
                residueA = chainA.residues[chainA.residue_indexes[i]]
                residueB = chainB.residues[chainB.residue_indexes[j]]
                print("Found", residueA, residueB, distance)
                interactions.append((residueA, residueB, distance))
                interaction_matrix[i,j] = 1.0
    return interactions, interaction_matrix

# CA distance calculation 
def get_CA_distance(residueA, residueB):
    return calculate_distance(residueA.atoms["CA"], residueB.atoms["CA"])

# Geting mean distance 
def get_mean_distance(residueA, residueB):
    return calculate_distance(residueA.get_mean_coordinate(), residueB.get_mean_coordinate())

# Creating distance matrixes 
def create_distance_matrix(chains, _chainA, _chainB, get_atom_distance):
    chainA = chains[_chainA]
    chainB = chains[_chainB]
    
    distance_matrix = np.zeros( (len(chainA.residues), len(chainB.residues)))
    for i, resnumA in enumerate(chainA.residue_indexes):
        residueA = chainA.residues[resnumA]
        for j, resnumB in enumerate(chainB.residue_indexes):
            residueB = chainB.residues[resnumB]
            distance = get_atom_distance(residueA, residueB)
            distance_matrix[i,j] = distance
    return distance_matrix

# Sub Matrix creation , Dina version
def create_fixedsize_submatrix(distmat_AB, sub_size, overlap):
    sub_mat = []
    rows, cols = distmat_AB.shape  
    for i in range(0, rows - sub_size +1, overlap):
        for j in range(0, cols - sub_size +1, overlap):
            #sub_matrix = distmat_AB[i:i+sub_size, j:j+sub_size]
            sub_matrix = Matrix(sub_size, sub_size)
            sub_matrix.matrix = distmat_AB[i:i+sub_size, j:j+sub_size]
            sub_matrix.i = i + 1  
            sub_matrix.j = j + 1  
            sub_mat.append(sub_matrix) # residue_indexes
    return sub_mat

def get_submatrix(distance_matrix,i,j, size):
    return distance_matrix[i:i+size,j:j+size]


        #sasa_pdb = result.write_pdb('2.pdb')       
def rsa(pdb_file, chains_CA): 
    structure = freesasa.Structure(pdb_file)
    result = freesasa.calc(structure)
    area_classes = freesasa.classifyResults(result, structure)
    for chain in chains_CA.values():
        if str(pdb_file).endswith(chain.chainID):
            for residue in chain.residues.values():
                print(chain.residue_indexes)
                all = result.residueAreas()               
                residue_sasa = all[chain.chainID][str(residue.resnum)].total
                residue.addRSA(residue_sasa)
                break  

c3_convert = {' ':0, 'S':0, 'T':0, 'H':1, 'G':1, 'I':1, 'E':2, 'B':2}


def dssp(chain_CA,pdb):
       # Sample coordinates
    
    for chain in chain_CA.values():
        length = len(chain.residues)
        all_atoms, all_coords = chain.get_all_atoms()
        atom = 5
        xyz = 4
    ## atoms should be 4 (N, CA, C, O) or 5 (N, CA, C, O, H)
        coord = torch.tensor(pydssp.read_pdbtext(open(pdb, 'r').read()))
    # hydrogene bond matrix
        hbond_matrix = pydssp.get_hbond_map(coord)
        print(hbond_matrix.shape) # should be (batch, length, length)
    # getting scondary struct 
        dssp_struct = pydssp.assign(coord, out_type='c3')
    ## output is batched np.ndarray of C3 annotation, like ['-', 'H', 'H', ..., 'E', '-']
    # To get secondary str. as index
        dssp_index = pydssp.assign(coord, out_type='index')
    ## 0: loop,  1: alpha-helix,  2: beta-strand
    # To get secondary str. as onehot representation
        dssp_onhot = pydssp.assign(coord, out_type='onehot')
    ## dim-0: loop,  dim-1: alpha-helix,  dim-2: beta-strand
        print(dssp_struct)
        print(hbond_matrix)


protein_data= {}
protein_sequences = []  
distance_matrices = []  
submatrices = []      
rsa_values = []  

def main():
    # Work directory
    work_dir = "/home/pc550/Documents/PPI_W_DLLM/workdir"
    processed_pdb_files = process_tgz_files_in_directory(work_dir) 
    chain_split_files = []
    # Looping over each pdb file in the directory 
    for pdb_file in processed_pdb_files:
        print("CA DISTANCE CALCULATION")
        print(pdb_file)
        chains_CA = parsePDB(pdb_file)
        [chainID1, chainID2] = chains_CA.keys()
        distance_matrix_CA__A_A = create_distance_matrix(chains_CA, chainID1, chainID1, get_CA_distance)
        distance_matrix_CA__B_B = create_distance_matrix(chains_CA, chainID2, chainID2, get_CA_distance)
        distance_matrix_CA__A_B = create_distance_matrix(chains_CA, chainID1, chainID2, get_CA_distance)
        interactions_CA__A_B, IM_CA__A_B = findInteractingResidues(chains_CA, chainID1, chainID2, distance_matrix_CA__A_B) # chains_CA)
        print(IM_CA__A_B)

        """ JS sub matrix creation
        [ row, col ] = IM_CA__A_B.shape
        
        for i in range(row-size):
            for j in range(col-size):
                interacting_pairs = np.sum(get_submatrix(IM_CA__A_B, i, j, size))
                if interacting_pairs > 1:
                    print(i,j, interacting_pairs)
                    print(get_submatrix(IM_CA__A_B, i, j, size))

        [ row, col ] = distance_matrix_CA__A_A.shape
        
        for i in range(row-size):
            for j in range(col-size):
                print(get_submatrix(distance_matrix_CA__A_A, i, j, size))
        """

        print("AA'S ATOMS DISTANCE CALCULATION")
        [chainID1, chainID2] = chains_CA.keys()
        distance_matrix_mean__A_A = create_distance_matrix(chains_CA, chainID1, chainID1, get_mean_distance)
        distance_matrix_mean__B_B = create_distance_matrix(chains_CA, chainID2, chainID2, get_mean_distance)
        distance_matrix_mean__A_B = create_distance_matrix(chains_CA, chainID1, chainID2, get_mean_distance)
        interactions_mean__A_B, IM_mean__A_B = findInteractingResidues(chains_CA, chainID1, chainID2, distance_matrix_mean__A_B) # chains_CA)

        print("SPEARMANR")
        print(stats.spearmanr(distance_matrix_CA__A_A.flatten(), distance_matrix_mean__A_A.flatten()))
        print(stats.spearmanr(distance_matrix_CA__B_B.flatten(), distance_matrix_mean__B_B.flatten()))
        print(stats.spearmanr(distance_matrix_CA__A_B.flatten(), distance_matrix_mean__A_B.flatten()))
        print(np.mean(abs(distance_matrix_CA__A_B-distance_matrix_mean__A_B)))

        # splitting chains and calculating the rsa on them 
        
        dir_name = os.path.dirname(pdb_file)
        splitted_files = splitPDBbyChain(pdb_file, dir_name)
        dir_name = os.path.dirname(pdb_file)
        chain_split_files.extend([os.path.join(dir_name, f) for f in os.listdir(dir_name) if f.endswith('.pdb') and f != os.path.basename(pdb_file)])
        
        for s in chain_split_files:
            rsa(s, chains_CA)
            dssp(chains_CA, s)
            break
        break         
if __name__ == "__main__":
    main()

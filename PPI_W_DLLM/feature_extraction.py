from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import tarfile
import os
import pickle
import numpy as np
import pandas as pd
from scipy import stats

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

class Chain:
    def __init__(self, chainID):
        self.chainID = chainID
        self.residues = {}
        self.residue_indexes = []       

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
    
class Residue:
    def __init__(self, aa, resnum):
        self.aa = aa
        self.resnum = resnum
        self.atoms = {}
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
    
    
def parsePDB(pdb_file):
    chains = {}
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
            if resnum not in chains[chainID].residues:
                chains[chainID].addResidue( aa, resnum )
            chains[chainID].residues[resnum].addAtom( atom_name, x, y, z ) # .addCA(aa, resnum, x, y, z)    
    return chains
"""
# Calculating each AA mean 
def calculating_aa_mean(pdb):
    chains_mean_points = {}
    with open(pdb, 'r') as f:
            current_chain = None
            current_residue = None
            residue_coordinates = []
            for line in f:
                if line.startswith("ATOM"):
                    aa = line[17:20]
                    resnum = int(line[22:26])
                    chainID = line[21]
                    if current_chain != chainID or current_residue != resnum:
                        if residue_coordinates:  
                            mean_x = sum(coord[0] for coord in residue_coordinates) / len(residue_coordinates)
                            mean_y = sum(coord[1] for coord in residue_coordinates) / len(residue_coordinates)
                            mean_z = sum(coord[2] for coord in residue_coordinates) / len(residue_coordinates)
                            mean_point = [mean_x, mean_y, mean_z]
                            if current_chain not in chains_mean_points:
                                chains_mean_points[current_chain] = Chain(current_chain)
                            chains_mean_points[current_chain].addCA(aa, current_residue, *mean_point)
                            residue_coordinates = []
                        current_chain = chainID
                        current_residue = resnum
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    residue_coordinates.append([x, y, z])
            if residue_coordinates:
                mean_x = sum(coord[0] for coord in residue_coordinates) / len(residue_coordinates)
                mean_y = sum(coord[1] for coord in residue_coordinates) / len(residue_coordinates)
                mean_z = sum(coord[2] for coord in residue_coordinates) / len(residue_coordinates)
                mean_point = [mean_x, mean_y, mean_z]
                if current_chain not in chains_mean_points:
                    chains_mean_points[current_chain] = Chain(current_chain)
                chains_mean_points[current_chain].addCA(aa, resnum, *mean_point)

    return chains_mean_points
"""

def calculate_distance(p1, p2):
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    dz = p1[2] - p2[2]
    distsq = pow(dx, 2) + pow(dy, 2) + pow(dz, 2)
    distance = np.sqrt(distsq)
    return distance


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

def get_CA_distance(residueA, residueB):
    return calculate_distance(residueA.atoms["CA"], residueB.atoms["CA"])

def get_mean_distance(residueA, residueB):
    return calculate_distance(residueA.get_mean_coordinate(), residueB.get_mean_coordinate())

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

#Dina version
def create_fixedsize_submatrix(distmat_AB, sub_size, overlap):
    sub_mat = []
    rows, cols = distmat_AB.shape  
    for i in range(0, rows - sub_size +1, overlap):
        for j in range(0, cols - sub_size +1, overlap):
            sub_matrix = distmat_AB[i:i+sub_size, j:j+sub_size]
            sub_mat.append((sub_matrix, i, j))
    return sub_mat

    '''
    unique_residues = {residue for interaction in interactions for residue in interaction[:2]}
    residue_to_index = {residue: i for i, residue in enumerate(unique_residues)}
    num_residues = len(unique_residues)
    distance_matrix = np.full((num_residues, num_residues), np.nan)
    for residueA, residueB, distance in interactions:
        indexA, indexB = residue_to_index[residueA], residue_to_index[residueB]
        distance_matrix[indexA, indexB] = distance
        distance_matrix[indexB, indexA] = distance  # Symmetric
        
    return distance_matrix
    '''
def get_submatrix(distance_matrix,i,j, size):
    return distance_matrix[i:i+size,j:j+size]

"""
# getting bost library issue : 
/home/pc550/miniconda3/lib/python3.12/site-packages/Bio/PDB/PDBParser.py:395: PDBConstructionWarning: Ignoring unrecognized record 'END' at line 1085
  warnings.warn(
/home/pc550/miniconda3/lib/python3.12/site-packages/Bio/PDB/DSSP.py:250: UserWarning: mkdssp: error while loading shared libraries: libboost_filesystem.so.1.73.0: cannot open shared object file: No such file or directory

  warnings.warn(err)
Traceback (most recent call last):
  File "/home/pc550/Documents/PPI_W_DLLM/YangLabIntern/PPI_W_DLLM/feature_extraction.py", line 339, in <module>
    main()
  File "/home/pc550/Documents/PPI_W_DLLM/YangLabIntern/PPI_W_DLLM/feature_extraction.py", line 273, in main
    rsa_data = calculate_rsa(pdb_file)
               ^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/pc550/Documents/PPI_W_DLLM/YangLabIntern/PPI_W_DLLM/feature_extraction.py", line 243, in calculate_rsa
    dssp = DSSP(model, pdb_file)
           ^^^^^^^^^^^^^^^^^^^^^
  File "/home/pc550/miniconda3/lib/python3.12/site-packages/Bio/PDB/DSSP.py", line 429, in __init__
    dssp_dict, dssp_keys = dssp_dict_from_pdb_file(in_file, dssp)
                           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/pc550/miniconda3/lib/python3.12/site-packages/Bio/PDB/DSSP.py", line 252, in dssp_dict_from_pdb_file
    raise Exception("DSSP failed to produce an output")
Exception: DSSP failed to produce an output
"""

"""
def calculate_rsa(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]
    chains = []
    for chain in model:
        chains.append(chain)
    rsa_data = {}
    for chain in chains:
        dssp = DSSP(model, pdb_file)
        rsa_data[chain.id] = {}
        for residue in chain:
            res_id = residue.id[1]
            rsa = dssp[(chain.id, res_id)]['RSA'] 
            rsa_data[chain.id][res_id] = rsa
    return rsa_data
"""

def main():
    # Work directory
    work_dir = "/home/pc550/Documents/PPI_W_DLLM/workdir"

    # Reading in the coordinates of each .pdb
    processed_pdb_files = process_tgz_files_in_directory(work_dir)

    #JS
    data1 = []
    data2 = []
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

        """
        rsa_data = calculate_rsa(pdb_file)
        print(rsa_data)
        """
        #data1.append(chains_CA)
        #distance_matrix = create_distance_matrix(interactions_CA)
        #print(distance_matrix)
        

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


        print(stats.spearmanr(distance_matrix_CA__A_A.flatten(), distance_matrix_mean__A_A.flatten()))
        print(stats.spearmanr(distance_matrix_CA__B_B.flatten(), distance_matrix_mean__B_B.flatten()))
        print(stats.spearmanr(distance_matrix_CA__A_B.flatten(), distance_matrix_mean__A_B.flatten()))
        print(np.mean(abs(distance_matrix_CA__A_B-distance_matrix_mean__A_B)))

        size = 7 
        overlap = 1 

        submatrices = create_fixedsize_submatrix(distance_matrix_CA__A_B, size, overlap)

        for idx, submatrix in enumerate(submatrices):
            print(f"Submatrix {idx+1}:")
            print(submatrix)
            """
            # SPEARMAN NOT WORKING CAUSE NOT THE RIGHT SIZE EVEN IF ITS VECTORIZED
            reshaped_submatrices = [submatrix.flatten() for submatrix in submatrices]
            stats.spearmanr(np.concatenate(reshaped_submatrices), distance_matrix_mean__A_B.flatten())
            print(np.mean(submatrices-distance_matrix_CA__A_B))
            """
    
        """ 
        chains_mean = calculating_aa_mean(pdb_file)
        distance_matrix_mean = create_distance_matrix(chains_mean, chainID2, chainID2)
        interactions_mean = findInteractingResidues(chains_mean, chainID2, chainID2, distance_matrix_mean) # chains_CA)  findInteractingFragments(chains_mean)
        data2.append(chains_mean)
        """

        break

if __name__ == "__main__":
    main()

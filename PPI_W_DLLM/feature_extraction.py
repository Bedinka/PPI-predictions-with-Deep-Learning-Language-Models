from Bio.PDB import PDBParser
import tarfile
import os
import numpy as np
import pandas as pd 
from scipy import stats
import freesasa
import torch 
import pydssp
import pickle

amino_acids_dict = {
    'ALA': 1,
    'ARG': 2,
    'ASN': 3,
    'ASP': 4,
    'CYS': 5,
    'GLN': 6,
    'GLU': 7,
    'GLY': 8,
    'HIS': 9,
    'ILE': 10,
    'LEU': 11,
    'LYS': 12,
    'MET': 13,
    'PHE': 14,
    'PRO': 15,
    'SER': 16,
    'THR': 17,
    'TRP': 18,
    'TYR': 19,
    'VAL': 20
}
work_dir = "/home/dina/Documents/PPI_WDLLM/workdir"
cache_dir = "cache"

# Function to load cached data
def load_cache(filename):
    try:
        with open(filename, "rb") as file:
            return pickle.load(file)
    except FileNotFoundError:
        return None

def save_cache(filename, data):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "wb") as file:
        pickle.dump(data, file)

def check_cache_and_calculate(cache_filename, calculate_function, *args, **kwargs):
    data = load_cache(cache_filename)
    if data is None:
        data = calculate_function(*args, **kwargs)
        save_cache(cache_filename, data)
    return data

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

    return splitted_files 

class Chain:
    def __init__(self, chainID, sample_number):
        self.samplenum = sample_number
        self.chainID = chainID
        self.residues = {}
        self.residue_indexes = []
        self.aa_dict = []    
        self.prot_id = []
        self.dssp_struct = []
        self.dssp_index = []
        self.dssp_onehot = []
        self.distance_matrices_CA = []
        self.distance_matrices_mean = []
        self.sub_res_index = []
        self.sub_res_name =[]
        self.distance_matrices_CA_AB = []
        self.distance_matrices_mean_AB = []
        self.mean_submatrices = []
        self.ca_submatrices = []
        self.interactions_CA =[]
        self.interactions_m_CA =[]
        self.interactions_mean =[]
        self.interactions_m_mean =[]
          

    def addResidue(self, aa, resnum):
        aResidue = Residue( aa, resnum )
        self.residues[resnum] = aResidue
        self.residue_indexes.append(resnum)
        if aa in amino_acids_dict:
            self.aa_dict.append(amino_acids_dict[aa])
        pass
    
    def addCA(self, aa, resnum, x, y, z):
        aResidue = Residue( aa, resnum, x, y, z)
        self.residues[resnum] = aResidue
        self.residue_indexes.append(resnum)
        pass

    def addCAMatrix (self, chains_CA , dist_mat , size, overlap ):
        aDistMatrix = Matrix(size, size)
        aDistMatrix.cadistmatrix = dist_mat
        self.distance_matrices_CA_AB = aDistMatrix.cadistmatrix
        #self.mean_submatrices = aDistMatrix.submatrixes(chains_CA, dist_mat , size, overlap)
        #CAREFULL: ADDING SUB MATRIX WITHOUT CREATING A MATRIX OBJECT
        self.ca_submatrices = create_fixedsize_submatrix(dist_mat, size, overlap , self.aa_dict)
        pass

    def addMeanMatrix (self, chains_CA , dist_mat , size, overlap):
        aDistMatrix = Matrix(size, size)
        aDistMatrix.meandistmatrix = dist_mat
        self.distance_matrices_mean_AB = aDistMatrix.meandistmatrix
        #self.mean_submatrices = aDistMatrix.submatrixes(chains_CA, dist_mat , size, overlap)
        #CAREFULL: ADDING SUB MATRIX WITHOUT CREATING A MATRIX OBJECT
        self.mean_submatrices = create_fixedsize_submatrix(dist_mat, size, overlap , self.aa_dict)
        self.sub_res_index , self.sub_res_name = sub_residuses (self.residues, self.residue_indexes, size)
        #NEED THE SUB nemas and indexes here too!!
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
        self.cadistmatrix= []
        self.meandistmatrix= []
        self.submatrix = []
        self.row = row
        self.col = col
        self.i=0 
        self.j=0
        self.residues_a = []
        self.residues_b = []

    def submatrixes( self, chains_CA, dist_mat , size, overlap ): # not in use 
        [chainID1, chainID2] = chains_CA.keys()
        submatrix = create_fixedsize_submatrix(dist_mat, size, overlap) 
        self.residuses_a , self.residues_b = sub_residuses( chains_CA, chainID1, chainID2, submatrix)
        self.submatrix =submatrix
        #print("Created fixed sized matrixes and got the starting residues")
        return submatrix
       
# Parsing the PDB files 
def parsePDB(pdb_file, sample_counter):
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
                chains[chainID] = Chain(chainID, sample_counter)
                if chainID == 'A':
                    chains[chainID].prot_id = first_part
                    #print(chains[chainID].prot_id)
                    sample_counter += 1
                else :
                    chains[chainID].prot_id = second_part
                    sample_counter += 1
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
                #print("Found", residueA, residueB, distance)
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
def create_fixedsize_submatrix(distmat, sub_size, overlap, aa_dict):
    sub_mat = []
    rows, cols = distmat.shape  
    for i in range(0, rows - sub_size + 1, overlap):
        for j in range(0, cols - sub_size + 1, overlap):
            sub_matrix = distmat[i:i+sub_size, j:j+sub_size]
            ''' Adding the index to the matrix 
            sub_matrix = np.zeros((sub_size + 1, sub_size + 1))
            sub_matrix[1:, 0] = aa_dict[i:i+sub_size]
            sub_matrix[0, 1:] = aa_dict[j:j+sub_size]
            sub_matrix[1:, 1:] = distmat[i:i+sub_size, j:j+sub_size]
            '''
            sub_mat.append(sub_matrix)
    '''for i in range(0, rows - sub_size +1, overlap):
        for j in range(0, cols - sub_size +1, overlap):
            sub_matrix = distmat[i:i+sub_size, j:j+sub_size]
            #sub_matrix = Matrix(sub_size, sub_size)
            #sub_matrix.matrix = distmat_AB[i:i+sub_size, j:j+sub_size]
            #sub_matrix.i = i + 1  
            #sub_matrix.j = j + 1  
            aa_slice = aa_dict[i:i+sub_size]  
            sub_mat.append((sub_matrix, aa_slice)) # residue_indexes'''
    return sub_mat

def get_submatrix(distance_matrix,i,j, size):
    return distance_matrix[i:i+size,j:j+size]

def sub_residuses( residues_list, residue_indexes, size ):
    sub_residues = []
    sub_names = []
    for i in range(len(residue_indexes) - size + 1):
            for j in range(len(residue_indexes) - size + 1):
                subr = [residue_indexes[i:i + size], residue_indexes[j:j + size]]
                sub_residues.append(subr)
                subn = [
                    [residues_list[resnum].aa for resnum in residue_indexes[i:i+size]],
                    [residues_list[resnum].aa for resnum in residue_indexes[j:j+size]]
                ]
                sub_names.append(subn)
    return sub_residues , sub_names

    
def rsa(pdb_file, chains_CA, work_dir): 
    try:
        #print(f"run FreeSASA on : {pdb_file}")
        structure = freesasa.Structure(pdb_file)
        result = freesasa.calc(structure)
        area_classes = freesasa.classifyResults(result, structure)
        result.write_pdb(os.path.basename(pdb_file))
        for chain in chains_CA.values():
            if str(pdb_file).endswith(chain.chainID):
                for residue in chain.residues.values():
                    all = result.residueAreas()               
                    residue_sasa = all[chain.chainID][str(residue.resnum)].total
                    residue.addRSA(residue_sasa)        
    except Exception as e:
        print(f"An error occurred while processing {pdb_file}: {e}")         

def dssp(chain_CA, pdb):
    for chain in chain_CA.values():
    ## atoms should be 4 (N, CA, C, O) or 5 (N, CA, C, O, H)
        coord = torch.tensor(pydssp.read_pdbtext(open(pdb, 'r').read()))
    # hydrogene bond matrix
        hbond_matrix = pydssp.get_hbond_map(coord)
        # should be (batch, length, length)
    # getting scondary struct 
        dssp_struct = pydssp.assign(coord, out_type='c3')
    # To get secondary str. as index
        dssp_index = pydssp.assign(coord, out_type='index')
        chain.dssp_indes = dssp_index
        ## 0: loop,  1: alpha-helix,  2: beta-strand
    # To get secondary str. as onehot representation
        dssp_onhot = pydssp.assign(coord, out_type='onehot')
        chain.dssp_onehot = dssp_onhot
        ## dim-0: loop,  dim-1: alpha-helix,  dim-2: beta-strand
        output_file =  work_dir + "/dssp/" + str(chain.prot_id) + ".dssp.txt"
        with open(output_file, 'w') as f:
            f.write(f"{hbond_matrix}\n{dssp_struct}\n{dssp_index}\n{dssp_onhot}\n")

def ca_dist_calc(chains_CA, size, overlap ):
    #print("CA DISTANCE CALCULATION") 
    [chainID1, chainID2] = chains_CA.keys()
    distance_matrix_CA__A_A = create_distance_matrix(chains_CA, chainID1, chainID1, get_CA_distance)
    chains_CA[chainID1].distance_matrices_CA = distance_matrix_CA__A_A
    distance_matrix_CA__B_B = create_distance_matrix(chains_CA, chainID2, chainID2, get_CA_distance)
    chains_CA[chainID2].distance_matrices_CA = distance_matrix_CA__B_B
    distance_matrix_CA__A_B = create_distance_matrix(chains_CA, chainID1, chainID2, get_CA_distance)
    chains_CA[chainID1].addCAMatrix(chains_CA, distance_matrix_CA__A_A, size , overlap ) # AB was chaged to AA
    chains_CA[chainID2].addCAMatrix(chains_CA, distance_matrix_CA__B_B, size , overlap )# AB was chaged to BB
    return distance_matrix_CA__A_B

def mean_dist_calc(chains_CA, size, overlap):
    #print("AA'S ATOMS DISTANCE CALCULATION")
    [chainID1, chainID2] = chains_CA.keys()
    distance_matrix_mean__A_A = create_distance_matrix(chains_CA, chainID1, chainID1, get_mean_distance)
    chains_CA[chainID1].distance_matrices_mean = distance_matrix_mean__A_A
    distance_matrix_mean__B_B = create_distance_matrix(chains_CA, chainID2, chainID2, get_mean_distance)
    chains_CA[chainID2].distance_matrices_mean = distance_matrix_mean__B_B
    distance_matrix_mean__A_B = create_distance_matrix(chains_CA, chainID1, chainID2, get_mean_distance)
    chains_CA[chainID1].addMeanMatrix(chains_CA, distance_matrix_mean__A_A, size , overlap )# AB was chaged to AA
    chains_CA[chainID2].addMeanMatrix(chains_CA, distance_matrix_mean__B_B, size , overlap )# AB was chaged to BB
    return distance_matrix_mean__A_B
    
def interacting_res (chains_CA, distance_matrix_CA__A_B , distance_matrix_mean__A_B):
    [chainID1, chainID2] = chains_CA.keys()
    interactions_CA__A_B, IM_CA__A_B = findInteractingResidues(chains_CA, chainID1, chainID2, distance_matrix_CA__A_B) 
    chains_CA[chainID1].interactions_CA = chains_CA[chainID2].interactions_CA = interactions_CA__A_B
    chains_CA[chainID1].interactions_m_CA = chains_CA[chainID2].interactions_m_CA = IM_CA__A_B

    interactions_mean__A_B, IM_mean__A_B = findInteractingResidues(chains_CA, chainID1, chainID2, distance_matrix_mean__A_B) 
    chains_CA[chainID1].interactions_mean = chains_CA[chainID2].interactions_mean = interactions_mean__A_B
    chains_CA[chainID1].interactions_m_mean = chains_CA[chainID2].interactions_m_mean = IM_mean__A_B
    return interactions_mean__A_B, IM_mean__A_B
"""
def spearman():
        print("SPEARMANR")
        print(stats.spearmanr(distance_matrix_CA__A_A.flatten(), distance_matrix_mean__A_A.flatten()))
        print(stats.spearmanr(distance_matrix_CA__B_B.flatten(), distance_matrix_mean__B_B.flatten()))
        print(stats.spearmanr(distance_matrix_CA__A_B.flatten(), distance_matrix_mean__A_B.flatten()))
        print(np.mean(abs(distance_matrix_CA__A_B-distance_matrix_mean__A_B)))   
"""
sample_counter = 1
def main(processed_sample, size):

    processed_pdb_files = [os.path.join(f"{work_dir}/interactions_001", file) for file in os.listdir(f"{work_dir}/interactions_001") if file.endswith('.pdb')]
   #processed_pdb_files = process_tgz_files_in_directory(work_dir)

    interacting_proteins = []
    overlap = 1 
    i = 1 # number of processed sample 

    for pdb_file in processed_pdb_files:
        chains_CA = parsePDB(pdb_file, sample_counter)
        [chainID1, chainID2] = chains_CA.keys()
        ca_dist = ca_dist_calc(chains_CA, size, overlap)
        mean_dist = mean_dist_calc(chains_CA, size, overlap)
        int_res = interacting_res(chains_CA, ca_dist, mean_dist )
        

        # splitting chains and calculating the rsa on them 
        dir_name = os.path.dirname(pdb_file)
        
        splitted_files = splitPDBbyChain(pdb_file, dir_name)

        for prot in splitted_files:
            rsa(prot, chains_CA, work_dir)
            dssp(chains_CA, prot)
        
        for chain in chains_CA.values():
            interacting_proteins.append(chain)      
        
        """cache_filename = os.path.join(cache_dir, f"{pdb_file}.cache")
        chains_CA = check_cache_and_calculate(cache_filename, parsePDB, pdb_file, sample_counter)

        cache_filename_ca = os.path.join(cache_dir, f"{pdb_file}_ca_dist.cache")
        cache_filename_mean = os.path.join(cache_dir, f"{pdb_file}_mean_dist.cache")
        cache_filename_int_res = os.path.join(cache_dir, f"{pdb_file}_int_res.cache")

        ca_dist = check_cache_and_calculate(cache_filename_ca, ca_dist_calc, chains_CA, size, overlap)
        mean_dist = check_cache_and_calculate(cache_filename_mean, mean_dist_calc, chains_CA, size, overlap)
        int_res = check_cache_and_calculate(cache_filename_int_res, interacting_res, chains_CA, ca_dist, mean_dist)

        dir_name = os.path.dirname(pdb_file)
        cache_filename_split = os.path.join(cache_dir, f"{pdb_file}_split.cache")
        splitted_files = check_cache_and_calculate(cache_filename_split, splitPDBbyChain, pdb_file, dir_name)

        for prot in splitted_files:
            cache_filename_rsa = os.path.join(cache_dir, f"{prot}_rsa.cache")
            rsa_data = check_cache_and_calculate(cache_filename_rsa, rsa, prot, chains_CA)
            cache_filename_dssp = os.path.join(cache_dir, f"{prot}_dssp.cache")
            dssp_data = check_cache_and_calculate(cache_filename_dssp, dssp, chains_CA, prot)"""

        for chain in chains_CA.values():
            interacting_proteins.append(chain)

        if i == processed_sample:
            break 
        
        #print(i)
        i += 1


    return interacting_proteins

if __name__ == "__main__":
    
    main()

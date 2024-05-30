#FEATURE EXTRACTION 
#Processing .pbd:
#Creating chain objects, creating CA and Mean distance matrix, creating CA and Mean submatrixes, .pickle from submatrixes
#Calculating RSA and DSSP 
#Creating .tsv from result with + and - samples

import tarfile
import os
import numpy as np
import pandas as pd 
from scipy import stats
import freesasa
import torch 
import pydssp


work_dir = "/home/dina/Documents/PPI_WDLLM/workdir"

class Chain:
    def __init__(self, chainID, sample_number):
        self.samplenum = sample_number
        self.chainID = chainID
        self.residues = {}
        self.residue_indexes = []
        self.aa_dict = []
        self.aa = []    
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
        self.int_prots = []
        self.interact = 0
          

    def addResidue(self, aa, resnum):
        aResidue = Residue( aa, resnum )
        self.residues[resnum] = aResidue
        self.residue_indexes.append(resnum)
        self.aa_dict.append(aa)
        pass
    
    def addCA(self, aa, resnum, x, y, z):
        aResidue = Residue( aa, resnum, x, y, z)
        self.residues[resnum] = aResidue
        self.residue_indexes.append(resnum)
        pass

    def addCAMatrix (self, chains_CA , dist_mat , size, overlap ):
        """aDistMatrix = Matrix(size, size)
        aDistMatrix.cadistmatrix = dist_mat"""
        self.distance_matrices_CA_AB = dist_mat
        #self.mean_submatrices = aDistMatrix.submatrixes(chains_CA, dist_mat , size, overlap)
        #CAREFULL: ADDING SUB MATRIX WITHOUT CREATING A MATRIX OBJECT
        self.ca_submatrices = create_fixedsize_submatrix(dist_mat, size, overlap )
        pass

    def addMeanMatrix (self, chains_CA , dist_mat , size, overlap):
        """aDistMatrix = Matrix(size, size)
        aDistMatrix.meandistmatrix = dist_mat"""
        self.distance_matrices_mean_AB = dist_mat
        #self.mean_submatrices = aDistMatrix.submatrixes(chains_CA, dist_mat , size, overlap)
        #CAREFULL: ADDING SUB MATRIX WITHOUT CREATING A MATRIX OBJECT
        self.mean_submatrices = create_fixedsize_submatrix(dist_mat, size, overlap )
        #self.sub_res_index , self.sub_res_name = sub_residuses (self.residues, self.residue_indexes, size)
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
def parsePDB(pdb_file_A, pdb_file_B, sample_counter):
    from Bio.PDB import PDBParser 
    from Bio.SeqUtils import seq1
    import warnings
    warnings.filterwarnings("ignore")
    chains = {}

    def extract_protein_id(pdb_file):
        with open(pdb_file, 'r') as f:
            title_count = 0
            for line in f:
                if line.startswith("TITLE"):
                    parts = line.split("(")
                    if len(parts) >= 2:
                        protein_id = parts[1].split(")")[0]
                        return protein_id.strip()
    def process_pdb_file(pdb_file, chain_label, sample_counter):
        parser = PDBParser()
        structure = parser.get_structure('structure', pdb_file)
        prot_id = extract_protein_id(pdb_file)
        local_counter = sample_counter  # Create a local copy of sample_counter
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM") == False:
                    continue
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                aa = line[17:20] 
                atom_name = line[12:16].strip()
                resnum = int(line[22:26]) 
                chainID = chain_label
                if chainID not in chains:
                    chains[chainID] = Chain(chainID, sample_counter)
                    chains[chainID].prot_id = prot_id
                    sample_counter += 1
                if resnum not in chains[chainID].residues:
                    chains[chainID].addResidue(aa, resnum)
                chains[chainID].residues[resnum].addAtom(atom_name, x, y, z)
            for model in structure:
                for chain in model:
                    chainID = chain.get_id()
                    if chainID == chain_label:
                        if chain_label not in chains:
                            chains[chain_label] = {'prot_id': prot_id, 'residues': {}, 'aa': []}
                        sequence = ''
                        for residue in chain:
                            if residue.get_id()[0] == ' ':
                                aa = seq1(residue.get_resname())
                                sequence += aa
                        chains[chain_label].aa.append(sequence)



    process_pdb_file(pdb_file_A, 'A', sample_counter)
    process_pdb_file(pdb_file_B, 'B', sample_counter)

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
def create_fixedsize_submatrix(distmat, sub_size, overlap):
    sub_mat = []
    rows, cols = distmat.shape  
    for i in range(0, rows - sub_size + 1, overlap):
        for j in range(0, cols - sub_size + 1, overlap):
            sub_matrix = distmat[i:i+sub_size, j:j+sub_size]
            sub_mat.append(sub_matrix)
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
        chain.dssp_struct = dssp_struct
    # To get secondary str. as index
        dssp_index = pydssp.assign(coord, out_type='index')
        chain.dssp_index = dssp_index
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

    
def interacting_res (chains_CA, distance_matrix_CA__A_B , distance_matrix_mean__A_B):
    [chainID1, chainID2] = chains_CA.keys()
    interactions_CA__A_B, IM_CA__A_B = findInteractingResidues(chains_CA, chainID1, chainID2, distance_matrix_CA__A_B) 
    chains_CA[chainID1].interactions_CA = chains_CA[chainID2].interactions_CA = interactions_CA__A_B
    chains_CA[chainID1].interactions_m_CA = chains_CA[chainID2].interactions_m_CA = IM_CA__A_B

    interactions_mean__A_B, IM_mean__A_B = findInteractingResidues(chains_CA, chainID1, chainID2, distance_matrix_mean__A_B) 
    chains_CA[chainID1].interactions_mean = chains_CA[chainID2].interactions_mean = interactions_mean__A_B
    chains_CA[chainID1].interactions_m_mean = chains_CA[chainID2].interactions_m_mean = IM_mean__A_B
    return interactions_mean__A_B, IM_mean__A_B

def create_positive_data_tsv(chains_CA, tsv_path):
    data = []
    [chainID1, chainID2] = chains_CA.keys()
    prot_id =  f"{chains_CA[chainID1].prot_id}_{chains_CA[chainID2].prot_id}" 
    aa = f"{chains_CA[chainID1].aa}_{chains_CA[chainID2].aa}" 
    dssp_struct = f"{chains_CA[chainID1].dssp_struct}_{chains_CA[chainID2].dssp_struct}" 
    dssp_index = f"{chains_CA[chainID1].dssp_index}_{chains_CA[chainID2].dssp_index}" 
    if chains_CA[chainID1].residue_indexes[-1] + 1 == chains_CA[chainID2].residue_indexes[0]:
        interact = 1
    else:
        interact = 0 
    data.append({
        'Protein ID': prot_id,
        'Interact': interact,
        'Residues' : aa,
        'DSSP Structure': dssp_struct,
        'DSSP Index': dssp_index

    })
    
    df = pd.DataFrame(data)
    if not os.path.isfile(tsv_path):
        df.to_csv(tsv_path, sep='\t', index=False)
    else:
        df.to_csv(tsv_path, sep='\t', index=False, header=False, mode='a')

def create_negative_data_tsv(interacting_proteins, tsv_path):
    data = []
    import random
    chains_CA = random.sample(interacting_proteins, 2)
    chainID1, chainID2 = chains_CA[0].chainID, chains_CA[1].chainID
    prot_id = f"{chains_CA[0].prot_id}_{chains_CA[1].prot_id}" 
    aa = f"{chains_CA[0].aa}_{chains_CA[1].aa}" 
    dssp_struct = f"{chains_CA[0].dssp_struct}_{chains_CA[1].dssp_struct}" 
    dssp_index = f"{chains_CA[0].dssp_index}_{chains_CA[1].dssp_index}" 
    if isinstance(chains_CA[0].int_prots, list) and chains_CA[1].prot_id in chains_CA[0].int_prots:
        interact = 1
    else:
        interact = 0 
        data.append({
            'Protein ID': prot_id,
            'Interact': interact,
            'Residues' : aa,
            'DSSP Structure': dssp_struct,
            'DSSP Index': dssp_index
        })
        
        df = pd.DataFrame(data)
        if not os.path.isfile(tsv_path):
                df.to_csv(tsv_path, sep='\t', index=False)
        else:
                df.to_csv(tsv_path, sep='\t', index=False, header=False, mode='a')

def matrix_pickle(chain_CA, pickle_ca_path, pickle_mean_path):
    import pickle
    for i in chain_CA:
        file_ca = 'm_ca_%s.pickle' % (chain_CA[i].prot_id)
        ca = os.path.join(pickle_ca_path, file_ca)

        if not os.path.exists(ca):
            with open(ca, 'wb') as f:
                pickle.dump(chain_CA[i].ca_submatrices, f)
            print(f"Created CA pickle file for {chain_CA[i].prot_id}")
        else:
            print(f"CA pickle file for {chain_CA[i].prot_id} already exists. Skipping.")

        file_mean = 'm_mean_%s.pickle' % (chain_CA[i].prot_id)
        mean = os.path.join(pickle_mean_path, file_mean)

        if not os.path.exists(mean):
            with open(mean, 'wb') as k:
                pickle.dump(chain_CA[i].mean_submatrices, k)
            print(f"Created mean pickle file for {chain_CA[i].prot_id}")
        else:
            print(f"Mean pickle file for {chain_CA[i].prot_id} already exists. Skipping.")

def find_pdb_file(protein_id, pdb_dir):
    dir = os.path.join(work_dir, pdb_dir)
    for file in os.listdir(dir):
        if protein_id in file:
            return os.path.join(dir, file)
    return None

sample_counter = 1


def main(processed_sample, size, tsv_path, pickle_ca_path, pickle_mean_path, pdb, positive_pairs_txt):
    print('Running Feature Extraction...')
    i = 1 # number of processed sample 
    interacting_proteins = []
    with open(positive_pairs_txt, 'r') as f:
        for line in f:
            pair = line.strip().split('\t')
            pdb1_protein_id = pair[0]
            pdb2_protein_id = pair[1]
            pdb1_file = find_pdb_file(pdb1_protein_id, pdb)
            pdb2_file = find_pdb_file(pdb2_protein_id, pdb)
            chains_CA = parsePDB(pdb1_file,pdb2_file, sample_counter) 
            #print(pdb_files)
            overlap = 1 
            print(pdb1_file, pdb2_file)
            [chainID1, chainID2] = chains_CA.keys()
            if chainID2 not in  chains_CA[chainID1].int_prots:
                chains_CA[chainID1].interact = 1
                chains_CA[chainID1].int_prots.append(chains_CA[chainID2].prot_id)
            if chainID1 not in chains_CA[chainID2].int_prots  :
                chains_CA[chainID2].interact = 1
                chains_CA[chainID2].int_prots.append(chains_CA[chainID1].prot_id)
            ca_dist_calc(chains_CA, size, overlap)
            mean_dist_calc(chains_CA, size, overlap)
            #int_res = interacting_res(chains_CA, ca_dist, mean_dist )
            
            import traceback
            try:
                #rsa(pdb_file, chains_CA, work_dir)
                dssp(chains_CA, pdb1_file)
                dssp(chains_CA, pdb2_file)
            except Exception as e:
                print(f"Error processing file : {e}")
                print(traceback.format_exc())
            
            for chain in chains_CA.values():
                if chain not in interacting_proteins:
                    interacting_proteins.append(chain)
            
            create_positive_data_tsv(chains_CA, tsv_path)

            matrix_pickle(chains_CA, pickle_ca_path, pickle_mean_path)

            if i == processed_sample:
                break 
            i += 1

        create_negative_data_tsv(interacting_proteins, tsv_path)
        print('Featur extraction : Done...')
    

if __name__ == "__main__":
    processed_sample = 4000
    size= 7
    tsv_path = 'bert_train_HuRI.tsv'
    pickle_ca_path= '/home/dina/Documents/PPI_WDLLM/Matrices_CA/new/'
    pickle_mean_path = '/home/dina/Documents/PPI_WDLLM/Matrices_Mean/new/'
    pdb = '/home/dina/Documents/PPI_WDLLM/workdir/Alpha_pdb'
    positive_pairs_txt = 'positive_pairs.txt'
    main(processed_sample, size, tsv_path, pickle_ca_path, pickle_mean_path, pdb,positive_pairs_txt)

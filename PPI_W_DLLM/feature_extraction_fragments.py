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
import traceback



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

    chain_files = []
    for chainID in chains:
        chain_filename = os.path.join(output_subdir, f"{basename}-{chainID}.pdb")
        with open(chain_filename, "w") as fout:
            fout.writelines(first_three_lines)
            for line in chains[chainID]:
                fout.write(line)
        chain_files.append(chain_filename)
    
    return chain_files[0], chain_files[1]

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
        self.rsa = []
        self.surface_residues = []
        self.surface_rsa = []
          

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
        #print(f"Added one rsa to {self.resnum}")
        #print(self.rsa_value)
         

# Parsing the PDB files 
def parsePDB(pdb_file, sample_counter):
    from Bio.PDB import PDBParser 
    from Bio.SeqUtils import seq1
    import warnings
    warnings.filterwarnings("ignore")
    chains = {}
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_file)
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
        for model in structure:
            for chain in model:
                chainID = chain.get_id()
                if chainID == 'A':
                    prot_id = first_part
                elif chainID == 'B':
                    prot_id = second_part
                else:
                    continue
                
                if chainID not in chains:
                    chains[chainID] = {'prot_id': prot_id, 'residues': {}, 'aa': []}
                    
                sequence = ''
                for residue in chain:
                    if residue.get_id()[0] == ' ':
                        aa = seq1(residue.get_resname())
                        sequence += aa
                
                chains[chainID].aa.append(sequence)
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

def get_atom_distance(residueA, residueB):
    # Assuming each residue object has a method to get the coordinates of the CA atom
    coordA = residueA.get_CA_coordinates()
    coordB = residueB.get_CA_coordinates()
    
    # Calculate the Euclidean distance between the two coordinates
    distance = np.linalg.norm(np.array(coordA) - np.array(coordB))
    return distance

# Creating distance matrixes 
def create_surface_distance_matrix(chainA, chainB, get_atom_distance):
    # Filter surface residues
    surface_residues_A = chainA.surface_residues
    surface_residues_B = chainB.surface_residues
    if surface_residues_A is None or surface_residues_B is None:
        print("One or both chains have no surface residues. Cannot create distance matrix.")
        return None
    distance_matrix = np.zeros((len(surface_residues_A), len(surface_residues_B)))
    for i, residueA in enumerate(surface_residues_A):
        for j, residueB in enumerate(surface_residues_B):
            distance = get_atom_distance(residueA, residueB)
            distance_matrix[i, j] = distance

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

    
def rsa(pdb_file, chain):
    try:
        structure = freesasa.Structure(pdb_file)
        result = freesasa.calc(structure)
        rsa_values = []
        surface_residues = []
        surface_rsa=[]
        threshold = 0.25  # This is a percentage value (25%)
        
        # Iterate through residues in the chain
        for residue in chain.residues.values():
            chain_id = chain.chainID
            if chain_id in result.residueAreas() and str(residue.resnum) in result.residueAreas()[chain_id]:
                # Get the RSA value for the residue
                residue_sasa = result.residueAreas()[chain_id][str(residue.resnum)].total
                residue.addRSA(residue_sasa)
                rsa_values.append(residue_sasa)
                
                if residue_sasa > threshold * 100:  
                    surface_residues.append(residue)
                    surface_rsa.append(residue_sasa)
                else:
                    print(f"not surface: {residue.resnum}")
            else:
                print(f"Chain {chain_id} or residue {residue.resnum} not found in the structure.")
        chain.rsa = rsa_values
        chain.surface_residues = surface_residues
        chain.surface_rsa = surface_rsa

        print(f"RSA values: {rsa_values}, Surface residues: {[res.resnum for res in surface_residues]}")
    except Exception as e:
        print(f"An error occurred: {e}")

def dssp(pdb_file, chain):
    try:
        with open(pdb_file, 'r') as pdb:
            coord = torch.tensor(pydssp.read_pdbtext(pdb.read()))

        hbond_matrix = pydssp.get_hbond_map(coord)
        dssp_struct = pydssp.assign(coord, out_type='c3')
        dssp_index = pydssp.assign(coord, out_type='index')
        dssp_onehot = pydssp.assign(coord, out_type='onehot')
        
        
        print(f"DSSP struct length: {len(dssp_struct)}, DSSP index length: {len(dssp_index)}")
        
        surface_residue_indices = [res.resnum - 1 for res in chain.surface_residues]
        print(f"Surface residue indices (0-based): {surface_residue_indices}")

        # Ensure indices are within bounds
        valid_surface_residue_indices = [idx for idx in surface_residue_indices if 0 <= idx < len(dssp_struct)]
        print(f"Valid surface residue indices: {valid_surface_residue_indices}")

        # Check for invalid indices
        invalid_indices = [idx for idx in surface_residue_indices if idx < 0 or idx >= len(dssp_struct)]
        if invalid_indices:
            print(f"Invalid surface residue indices: {invalid_indices}")

        surface_dssp_struct = [dssp_struct[idx] for idx in valid_surface_residue_indices]
        surface_dssp_index = [dssp_index[idx] for idx in valid_surface_residue_indices]
        surface_dssp_onehot = dssp_onehot[valid_surface_residue_indices]
        
        chain.dssp_struct = surface_dssp_struct
        chain.dssp_index = surface_dssp_index
        chain.dssp_onehot = surface_dssp_onehot
        
        # Detailed log of DSSP data for surface residues
        for idx in valid_surface_residue_indices:
            print(f"Index: {idx}, DSSP struct: {dssp_struct[idx]}, DSSP index: {dssp_index[idx]}")
        
    except IndexError as e:
        print(f"Index error in DSSP processing: {e}")
        print(traceback.format_exc())
    except Exception as e:
        print(f"An unexpected error occurred in DSSP processing: {e}")
        print(traceback.format_exc())


def ca_dist_calc(chains_CA, size, overlap ):
    #print("CA DISTANCE CALCULATION") 
    [chainID1, chainID2] = chains_CA.keys()
    distance_matrix_CA__A_A = create_surface_distance_matrix(chains_CA[chainID1], chains_CA[chainID1], get_CA_distance)
    chains_CA[chainID1].distance_matrices_CA = distance_matrix_CA__A_A
    distance_matrix_CA__B_B = create_surface_distance_matrix(chains_CA[chainID2], chains_CA[chainID2], get_CA_distance)
    chains_CA[chainID2].distance_matrices_CA = distance_matrix_CA__B_B
    #distance_matrix_CA__A_B = create_surface_distance_matrix(chains_CA[chainID1], chains_CA[chainID2], get_CA_distance)
    chains_CA[chainID1].addCAMatrix(chains_CA, distance_matrix_CA__A_A, size , overlap ) # AB was chaged to AA
    chains_CA[chainID2].addCAMatrix(chains_CA, distance_matrix_CA__B_B, size , overlap )# AB was chaged to BBs


def mean_dist_calc(chains_CA, size, overlap):
    #print("AA'S ATOMS DISTANCE CALCULATION")
    [chainID1, chainID2] = chains_CA.keys()
    distance_matrix_mean__A_A = create_surface_distance_matrix(chains_CA, chainID1, chainID1, get_mean_distance)
    chains_CA[chainID1].distance_matrices_mean = distance_matrix_mean__A_A
    distance_matrix_mean__B_B = create_surface_distance_matrix(chains_CA, chainID2, chainID2, get_mean_distance)
    chains_CA[chainID2].distance_matrices_mean = distance_matrix_mean__B_B
    distance_matrix_mean__A_B = create_surface_distance_matrix(chains_CA, chainID1, chainID2, get_mean_distance)
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
"""
def spearman():
        print("SPEARMANR")
        print(stats.spearmanr(distance_matrix_CA__A_A.flatten(), distance_matrix_mean__A_A.flatten()))
        print(stats.spearmanr(distance_matrix_CA__B_B.flatten(), distance_matrix_mean__B_B.flatten()))
        print(stats.spearmanr(distance_matrix_CA__A_B.flatten(), distance_matrix_mean__A_B.flatten()))
        print(np.mean(abs(distance_matrix_CA__A_B-distance_matrix_mean__A_B)))   
"""
def create_positive_data_tsv(chains_CA, tsv_path):
    data = []
    [chainID1, chainID2] = chains_CA.keys()

    if chains_CA[chainID2].prot_id in chains_CA[chainID1].int_prots:
        interact = 1
        prot_id =  f"{chains_CA[chainID1].prot_id}_{chains_CA[chainID2].prot_id}" 
        aa = f"{[res.aa for res in chains_CA[chainID1].surface_residues]}_{[res.aa for res in chains_CA[chainID2].surface_residues]}" 
        dssp_struct = f"{chains_CA[chainID1].dssp_struct}_{chains_CA[chainID2].dssp_struct}" 
        dssp_index = f"{chains_CA[chainID1].dssp_index}_{chains_CA[chainID2].dssp_index}" 
        rsa = f"{chains_CA[chainID1].surface_rsa}_{chains_CA[chainID2].surface_rsa}" 

        data.append({
        'Protein ID': prot_id,
        'Interact': interact,
        'Residues' : aa,
        'DSSP Structure': dssp_struct,
        'DSSP Index': dssp_index,
        'RSA' : rsa

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
    if chains_CA[1].prot_id in chains_CA[0].int_prots:
        interact = 1
    else:
        interact = 0 
        prot_id =  f"{chains_CA[0].prot_id}_{chains_CA[1].prot_id}" 
        aa = f"{[res.aa for res in chains_CA[0].surface_residues]}_{[res.aa for res in chains_CA[1].surface_residues]}" 
        dssp_struct = f"{chains_CA[0].dssp_struct}_{chains_CA[1].dssp_struct}" 
        dssp_index = f"{chains_CA[0].dssp_index}_{chains_CA[1].dssp_index}" 
        rsa = f"{chains_CA[0].surface_rsa}_{chains_CA[1].surface_rsa}" 
        data.append({
            'Protein ID': prot_id,
            'Interact': interact,
            'Residues' : aa,
            'DSSP Structure': dssp_struct,
            'DSSP Index': dssp_index,
            'RSA' : rsa
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
        chain_CA[i].ca_submatrices = None
        chain_CA[i].mean_submatrices = None
        print(f"Deleted ca_submatrices and mean_submatrices for {chain_CA[i].prot_id}")

sample_counter = 1


def main(processed_sample, size, tsv_path, pickle_ca_path, pickle_mean_path, pdb):
    print('Running Feature Extraction...')
    processed_pdb_files = [os.path.join(pdb, file) for file in os.listdir(pdb) if file.endswith('.pdb')]
    import random
    random.shuffle(processed_pdb_files)
    interacting_proteins = []
    overlap = 1 
    i = 1 # number of processed sample 

    for pdb_file in processed_pdb_files:
        print(pdb_file)
        chains_CA = parsePDB(pdb_file, sample_counter)
        [chainID1, chainID2] = chains_CA.keys()
        if chainID2 not in  chains_CA[chainID1].int_prots:
            chains_CA[chainID1].interact = 1
            chains_CA[chainID1].int_prots.append(chains_CA[chainID2].prot_id)
        if chainID1 not in chains_CA[chainID2].int_prots  :
            chains_CA[chainID2].interact = 1
            chains_CA[chainID2].int_prots.append(chains_CA[chainID1].prot_id)

        # splitting chains and calculating the rsa on them 
        dir_name = os.path.dirname(pdb_file)
        
        prot1, prot2 = splitPDBbyChain(pdb_file, dir_name)
        import traceback
        try:
            rsa( prot1,chains_CA[chainID1])
            dssp(prot1,chains_CA[chainID1])
        except Exception as e:
            print(f"Error processing file {prot1}: {e}")
            print(traceback.format_exc())
        try:
            rsa( prot2,chains_CA[chainID2])
            dssp(prot2,chains_CA[chainID2])
        except Exception as e:
            print(f"Error processing file {prot2}: {e}")
            print(traceback.format_exc())

        ca_dist_calc(chains_CA, size, overlap)
        #mean_dist_calc(chains_CA, size, overlap)
        #int_res = interacting_res(chains_CA, ca_dist, mean_dist )

        
        for chain in chains_CA.values():
            if chain not in interacting_proteins:
                interacting_proteins.append(chain)
        
        create_positive_data_tsv(chains_CA, tsv_path)

        matrix_pickle(chains_CA, pickle_ca_path, pickle_mean_path)

        if i == processed_sample:
            break
        i += 1
    for k in range(processed_sample):
        create_negative_data_tsv(interacting_proteins, tsv_path)
    print('Featur extraction : Done...')
    

if __name__ == "__main__":
    processed_sample = 10
    size= 7
    tsv_path = '/Users/baledi/ProtSeq Project /fragment.tsv'
    pickle_ca_path= '/Users/baledi/ProtSeq Project /Matrices_CA'
    pickle_mean_path = '/Users/baledi/ProtSeq Project /Matrices_Mean/'
    pdb = '/Users/baledi/ProtSeq Project /interactions_001'
    main(processed_sample, size, tsv_path, pickle_ca_path, pickle_mean_path, pdb)

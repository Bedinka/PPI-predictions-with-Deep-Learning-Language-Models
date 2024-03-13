from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import tarfile
import os
import numpy as np
import pandas as pd
from scipy import stats
import freesasa

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

def splitPDBbyChain(filename, output_dir):
    chains = {}
    with open(filename, 'r') as f:
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
            for line in chains[chainID]:
                fout.write(line)
    

class Chain:
    def __init__(self, chainID):
        self.chainID = chainID
        self.residues = {}
        self.residue_indexes = []       

    def addResidue(self, aa, resnum, chainID):
        aResidue = Residue( aa, resnum, chainID )
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
    def __init__(self, aa, resnum , chaiID ):
        self.aa = aa
        self.resnum = resnum
        self.chain = chaiID
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
        

        """"
        residue_structure = freesasa.Structure(pdb_file)
        for atom_name, (x, y, z) in self.atoms.items():
            residue_structure.addAtom(atom_name, x, y, z)
        result = freesasa.calc(residue_structure)
        self.rsa_value = result.totalArea()
        return rsa
        """
    
    
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
                chains[chainID].addResidue( aa, resnum , chainID)
            chains[chainID].residues[resnum].addAtom( atom_name, x, y, z ) # .addCA(aa, resnum, x, y, z)    
    return chains

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
            #sub_matrix = distmat_AB[i:i+sub_size, j:j+sub_size]
            sub_matrix = Matrix(sub_size, sub_size)
            sub_matrix.matrix = distmat_AB[i:i+sub_size, j:j+sub_size]
            sub_matrix.i = i + 1  
            sub_matrix.j = j + 1  
            sub_mat.append(sub_matrix) # residue_indexes
    return sub_mat

def get_submatrix(distance_matrix,i,j, size):
    return distance_matrix[i:i+size,j:j+size]

def rsa(pdb_files, chains_CA): 
    for pdb_file in pdb_files:
        structure = freesasa.Structure(pdb_file)
        result = freesasa.calc(structure)
        area_classes = freesasa.classifyResults(result, structure)
        for chain in chains_CA.values():
            for residue in chain.residues.values():
                print(chain.residue_indexes)
                all = result.residueAreas()
                residue_sasa = all[residue.chain][str(residue.resnum)].total
                for x in all[residue.chain]:
                    print(all[residue.chain][x].total)
                residue.addRSA(residue_sasa)
        break
        #sasa_pdb = result.write_pdb('2.pdb')       
    

    """ 
    Added one rsa to 14
140.62440188087328
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64]
115.51167668357688
113.77586293720269
130.8781912231541
53.48222976061239
24.55479769595844
47.312985416460144
137.25484826705383
16.212368839951854
44.8651779151281
90.81730768696035
100.90717367016215
34.291503137962266
118.94719638011956
140.62440188087328
164.45513059610005
Added one rsa to 15
164.45513059610005
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64]
Traceback (most recent call last):
  File "/home/pc550/Documents/PPI_W_DLLM/YangLabIntern/PPI_W_DLLM/feature_extraction.py", line 370, in <module>
    main()
  File "/home/pc550/Documents/PPI_W_DLLM/YangLabIntern/PPI_W_DLLM/feature_extraction.py", line 364, in main
    rsa(chain_split_files, chains_CA)    
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/pc550/Documents/PPI_W_DLLM/YangLabIntern/PPI_W_DLLM/feature_extraction.py", line 285, in rsa
    residue_sasa = all[residue.chain][str(residue.resnum)].total
                   ~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^
    """
def main():
    # Work directory
    work_dir = "/home/pc550/Documents/PPI_W_DLLM/workdir"
    
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
        #residue = Residue(aa,resnum)
        #rsa_value = residue.calculate_RSA(pdb_file)   
        chain_split_files = []
        for pdb_file in processed_pdb_files:
            dir_name = os.path.dirname(pdb_file)
            splitPDBbyChain(pdb_file, dir_name)
            dir_name = os.path.dirname(pdb_file)
            chain_split_files.extend([os.path.join(dir_name, f) for f in os.listdir(dir_name) if f.endswith('.pdb') and f != os.path.basename(pdb_file)])
        
        rsa(chain_split_files, chains_CA)    
        
            
        break

if __name__ == "__main__":
    main()

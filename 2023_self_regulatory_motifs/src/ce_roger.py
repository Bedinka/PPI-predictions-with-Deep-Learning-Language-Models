#!/usr/bin/env python3

import random
import os
import shutil
import numpy as np
import subprocess

class CEAlignmentMatrix:
    def __init__(self, directory, n):
        # Initialize the CEAlignmentMatrix class.
        self.directory = directory
        self.prot_label = random.sample(os.listdir(directory), n)
        #self.prot_label = ['AF-Q40240-F1-model_v4_7-156_PF01620.19.pdb', 'AF-P22285-F1-model_v4_44-176_PF01620.19.pdb', 'AF-P56167-F1-model_v4_2-28_PF01620.19.pdb', 'AF-Q40962-F1-model_v4_8-148_PF01620.19.pdb', 'AF-Q40237-F1-model_v4_7-48_PF01620.19.pdb', 'AF-P56166-F1-model_v4_94-140_PF01620.19.pdb', 'AF-P56164-F1-model_v4_7-145_PF01620.19.pdb', 'AF-P43215-F1-model_v4_2-130_PF01620.19.pdb', 'AF-Q40237-F1-model_v4_44-175_PF01620.19.pdb', 'AF-P22286-F1-model_v4_7-169_PF01620.19.pdb']
        self.matrix = np.zeros((n, n))

    #3 4
    #AF-Q40960-F1-model_v4_7-39_PF01620.19.pdb AF-P56166-F1-model_v4_94-140_PF01620.19.pdb
    #4 7
    #AF-P56166-F1-model_v4_94-140_PF01620.19.pdb 
    # AF-Q40960-F1-model_v4_7-39_PF01620.19.pdb

    def fill_rmsd(self):
        n = len(self.prot_label)
        for i in range(n):
            for j in range(i+1,n):
                # print(i,j)
                # print(self.prot_label[i], self.prot_label[j])
                # Calculate RMSD values between protein structures and fill the matrix.
                rmsd = self.calculate_rmsd(self.directory+'/'+self.prot_label[i], self.directory+'/'+self.prot_label[j])
                print(rmsd)
                self.matrix[i, j] = rmsd
                self.matrix[j, i] = self.matrix[i, j]
        return

    def find_best_model(self):
        # Find the best model based on the RMSD matrix.
        l = np.sum(self.matrix, axis=1)
        print(l)
        min_idx = np.argmin(np.abs(l - np.min(l)))
        print(self.prot_label[min_idx], np.min(l))
        self.best_model = self.prot_label[min_idx]

    def print_matrix(self):
        for row in self.matrix:
            print(row)

    def calculate_rmsd(self, PDB1, PDB2):
        # Calculate RMSD between two protein structures.
        # PDB1 = "/home/roger/2023_self_regulatory_motifs/SwissProt/domainStructures/PF01620.19/AF-P56166-F1-model_v4_94-140_PF01620.19.pdb" 
        # PDB2 = "/home/roger/2023_self_regulatory_motifs/SwissProt/domainStructures/PF01620.19/AF-Q40960-F1-model_v4_7-39_PF01620.19.pdb"
        print("python3 ../src/calculate_rmsd.py %s %s" % (PDB1, PDB2))

        result = subprocess.run(["python3", "../src/calculate_rmsd.py", PDB1, PDB2], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # print(result)
        # print(result.stderr)
        if (result.returncode==0):
            print("ok")
            rmsd = result.stdout
        else:
            rmsd = 100
            
        return rmsd    


def process_directory(directory, n=30):
    print('-------------------')
    print(directory)
    best_file = os.path.join(directory, 'best_model.pdb')
    # If best_model already exists delete it
    if os.path.exists(best_file):
        os.remove(best_file)
    protein_files = os.listdir(directory)
    # If there is only one structure this will be the best model
    if len(protein_files) == 1:
        best_model = protein_files[0]
        print('n = 1')
    else:
        # CHECKING THAT THERE ARE ENOUGH STRUCTURES ON THE FOLDER
        if len(protein_files) < n:
            n = len(protein_files)
        print('n =', n)

        m = CEAlignmentMatrix(directory, n)
        print('Object created')
        m.fill_rmsd()
        print('Matrix filled')
        m.print_matrix()
        m.find_best_model()
        
        # Get the best model filename
        best_model = m.best_model
    print('model found')
    # Copy the best model file and rename it
    original_path = os.path.join(directory, best_model)
    # Rename the file to 'best_model.pdb'
    new_path = os.path.join(directory, 'best_model.pdb')
    
    shutil.copyfile(original_path, new_path)
    
    print('Best model saved:', new_path)


pdbdir = '/home/roger/2023_self_regulatory_motifs/SwissProt/domainStructures'
c=0
# Process directories one by one
for entry in os.scandir(pdbdir):
   if entry.is_dir():
# for i in range(1000):
        # entrypath ='/home/roger/2023_self_regulatory_motifs/SwissProt/domainStructures/PF01620.19'
        process_directory(entry.path, n=30)
        c+=1
        print('c:', c)

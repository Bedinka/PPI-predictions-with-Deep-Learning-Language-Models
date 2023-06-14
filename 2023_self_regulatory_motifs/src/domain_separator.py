#!/usr/bin/env python3

import shelve
import os
from Bio import PDB

# CustomSelect class for selective output
class CustomSelect(PDB.Select):
    def __init__(self, start_index, end_index):
        self.start_index = start_index
        self.end_index = end_index

    def accept_residue(self, residue):
        residue_index = residue.get_id()[1]
        if self.start_index <= residue_index <= self.end_index:
            return 1  # Accept the residue
        else:
            return 0  # Reject the residue

x = shelve.open('domdict.shelve')

parentDir = 'domainStructures'

c = 0
for p in x.keys():
    for h in x[p]:
        for hsp in h:
            prot_name = p
            domain_id = h.accession
            domain_start = hsp.query_range[0] + 1
            domain_end = hsp.query_range[1]
            # create the domain directory
            path = os.path.join(parentDir, domain_id)
            # creating new pdb file
            parser = PDB.PDBParser()
            structure = parser.get_structure(prot_name, 'data/' + prot_name + '.pdb')
            io = PDB.PDBIO()
            io.set_structure(structure)
            # Selective output
            select = CustomSelect(start_index=domain_start, end_index=domain_end)
            output_file = os.path.join(path, f'{prot_name}_{domain_start}-{domain_end}_{domain_id}.pdb')
            try:
                io.save(output_file, select=select)
                print(f"Saved protein {prot_name} in {output_file}")
            except FileNotFoundError:
                # create the directory if it is not found
                os.mkdir(path, 0o777)
                print(f"Created directory: {path}")
                io.save(output_file, select=select)
                print(f"Saved protein {prot_name} in {output_file}")

            c += 1
            print(c)
            print('-----------------------')

x.close()

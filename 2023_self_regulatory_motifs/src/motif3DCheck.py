#!/usr/bin/env python3

import os, csv
from Bio.PDB import CEAligner, PDBParser
from Bio.PDB.cealign import _RESID_SORTER

# GENERATE A DICT THAT CONTAINS DOMAINS AS KEYS AND MOTIFS AS VALUES
motifDict = dict()
with open('SwissProt_results.csv') as csvMotifs:
    readMotifs= csv.reader(csvMotifs)
    h = next(readMotifs)
    for row in readMotifs:
        if row[2] not in motifDict.keys():
            motifDict[row[2]] = [row]
        else:
            motifDict[row[2]].append(row)
    print('csv parsed')
print(h)
# for i in motifDict['PF02980.19']:
#     print(i)
#     print()
pdbdir = '/home/roger/2023_self_regulatory_motifs/SwissProt/domainStructures'
c = 0
# Process directories one by one
for entry in os.scandir(pdbdir):
    if entry.is_dir():
        # Load the reference structure
        reference_structure = PDBParser().get_structure("reference", os.path.join(entry.path, "best_model.pdb"))
        if entry.name not in motifDict.keys():
            continue
        # for loop that saves the positions of the motifs after alignment
        for idx, motif in enumerate(motifDict[entry.name]):
            print(motif)
            # AF-Q9ZW31-F1-model_v4_49-152_PF00011.24.pdb
            structure_entry= motif[1] +'_'+ motif[5] +'-'+ motif[6] +'_'+ motif[2] + '.pdb'
            print(structure_entry)
            # Load the structure to align
            structure = PDBParser().get_structure("structure", os.path.join(entry.path, structure_entry))
            
            # Align structure to the reference
            aligner1 = CEAligner()
            aligner1.set_reference(reference_structure)
            aligner1.align(structure)
            # Save the translation and rotation (code line is missing, add appropriate code here)
            moves = aligner1.get_guide_coord_from_structure(structure)
            # Align full structure to moved structure
            full_structure = PDBParser().get_structure("full_structure", os.path.join('/home/roger/2023_self_regulatory_motifs/SwissProt/data/', structure_entry.split('_')[0]+'_'+structure_entry.split('_')[1]+'.pdb'))
            aligner2 = CEAligner()
            aligner2.set_reference(structure)
            aligner2.align(full_structure)
            should_be_moves = aligner2.refcoord
            assert(should_be_moves==moves)
            motifStart = int(motif[7])
            motifEnd = int(motif[8])
            motifPosition = aligner2.get_guide_coord_from_structure(full_structure)[motifStart:motifEnd]
            motifDict[entry.name][idx].append(motifPosition)
            print('\n----------------------------')
        # for loop that loops over the domains and checks if the motids are clustered
        for motif in motifDict[entry.name]:
            pass
        c += 1
        print('c:', c)
        break

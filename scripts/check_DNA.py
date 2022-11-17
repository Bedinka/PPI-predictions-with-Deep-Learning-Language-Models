#!/usr/bin/env python
"""
Check for PDB files containing DNA
"""



print(len(pdb_id))
dna=0
rna=0
dna_rna=0
for pdb in pdb_id:
    f=open(f"pdb_new_files/{pdb}/{pdb}.pdb")
    for lines in f.readlines():
        if "MOLECULE: DNA" in lines:
            print(pdb, lines)
            dna+=1
        else:
            continue
            
    for lines in f.readlines():
        if "MOLECULE: RNA" in lines:
            print(pdb, lines)
            rna+=1
        else:
            continue
            
    for lines in f.readlines():
        if "MOLECULE: DNA" and "MOLECULE: RNA" in lines:
            print(pdb, lines)
            dna_rna+=1
        else:
            continue

print(f"total number of structures with DNA/RNA: {dna_rna}")
print(f"total number of structures with DNA: {dna}")
print(f"total number of structures with RNA: {rna}")

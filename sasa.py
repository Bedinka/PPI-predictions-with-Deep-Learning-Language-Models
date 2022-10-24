import pandas as pd 
import os 
import freesasa
import glob 


"""
This script will estimate the solvent accessible surface area at protein and residual levels 
"""

#Extract PDB IDs and PDB chains 
pdbs = []
pdb_dict = {}

for i in os.listdir():
    if "." in i:
        continue
    else:
        pdbs.append(i)
        
for i in pdbs: 
    sel = glob.glob(f"{i}/*.pdb")
    for a in range(len(sel)):
        pdb_dict[i] = [sel[a-1].split(".")[0], sel[a].split(".")[0]]

#Calculate SASA values 
monomer=[]

for key,value in pdb_dict.items():
    
    chain1_rname, chain1_rnumber, chain1_relative, chain1_total = [], [], [], []
    
    c1_delta_relative, c1_delta_total= [], []

    chain2_rname, chain2_rnumber, chain2_relative, chain2_total= [], [], [], []
    
    c2_delta_relative, c2_delta_total = [], []

    comp1_rname, comp1_rnumber, comp1_relative, comp1_total = [], [], [], []
    
    comp2_rname, comp2_rnumber, comp2_relative, comp2_total = [], [], [], []

    combination1, combination2 = [], []

    
    
    a=value[0].split("/")[1] #Select chain 1
    b=value[1].split("/")[1]  #select chain 2
    id_a = value[0].split("-")[1].split(".")[0]
    id_b = value[1].split("-")[1].split(".")[0]
    pdb=key
    
    #Estimate SASA for chain 1
    structure=freesasa.Structure(f"{pdb}/{a}.pdb") 
    cal=freesasa.calc(structure)  
    total=cal.totalArea() 
    residues=cal.residueAreas()  

    ra=residues[list(residues.keys())[0]]

    for i in ra:
        ra=residues[list(residues.keys())[0]][i]
        chain1_rnumber.append(ra.residueNumber)
        chain1_rname.append(ra.residueType)
        chain1_relative.append(ra.relativeTotal)
        chain1_total.append(ra.total)
        
    #Estimating SASA for chain 2
    
    structure2=freesasa.Structure(f"{pdb}/{b}.pdb") 
    cal2=freesasa.calc(structure2)  
    total=cal2.totalArea() 
    residues=cal2.residueAreas()  

    rb=residues[list(residues.keys())[0]]

    for i in rb:
        rb=residues[list(residues.keys())[0]][i]
        chain2_rnumber.append(rb.residueNumber)
        chain2_rname.append(rb.residueType)
        chain2_relative.append(rb.relativeTotal)
        chain2_total.append(rb.total)
    

    

    #Estimating SASA for complex structure
    os.system("cat {}/{}.pdb {}/{}.pdb > {}/complex.pdb".format(pdb,a,pdb,b,pdb))
    complex_structure=freesasa.Structure("{}/complex.pdb".format(pdb))
    complex_cal=freesasa.calc(complex_structure)
    total_complex=complex_cal.totalArea()
    residues_complex=complex_cal.residueAreas()
    
    if len(residues_complex.keys()) != 1: 
        c1=residues_complex[list(residues_complex.keys())[0]]  
        c2=residues_complex[list(residues_complex.keys())[1]]

        for x in c1:
            c1=residues_complex[list(residues_complex.keys())[0]][x]
            comp1_rnumber.append(c1.residueNumber)
            comp1_rname.append(c1.residueType)
            comp1_relative.append(c1.relativeTotal)
            comp1_total.append(c1.total)

        for x in c2:
            c2=residues_complex[list(residues_complex.keys())[1]][x]
            comp2_rnumber.append(c2.residueNumber)
            comp2_rname.append(c2.residueType)
            comp2_relative.append(c2.relativeTotal)
            comp2_total.append(c2.total)


        for num in range(len(chain1_total)-1):
            c1_delta_total.append(chain1_total[num] - comp1_total[num])
            c1_delta_relative.append(chain1_relative[num] - comp1_relative[num])

        for i in range(len(chain2_total)-1):
            c2_delta_total.append(chain2_total[i] - comp2_total[i])
            c2_delta_relative.append(chain2_relative[i] - comp2_relative[i])

        #Saving data to csv files in seperate directories 

        combination1=list(zip(chain1_rname,chain1_rnumber,chain1_relative,chain1_total,
                        comp1_relative,comp1_total,c1_delta_relative,c1_delta_total))

        combination2=list(zip(chain2_rname,chain2_rnumber,chain2_relative,chain2_total,
                        comp2_relative,comp2_total,c2_delta_relative,c2_delta_total))



        pd.DataFrame(combination1,columns=["Residue name","Residue number","M.relative SASA","M.absolute SASA",
                                     "C.relative SASA","C.absolute SASA","Delta Relative","Delta Absolute"]).to_csv(f"{pdb}/SASA_{id_a}.csv",header=True)

        pd.DataFrame(combination2,columns=["Residue name","Residue number","M.relative SASA","M.absolute SASA",
                                     "C.relative SASA","C.absolute SASA","Delta Relative","Delta Absolute"]).to_csv(f"{pdb}/SASA_{id_b}.csv",header=True)

    else:
        monomer.append(pdb)



#!/usr/bin/env python
import glob
import os 
import freesasa
import pandas as pd
import numpy as np 



"""
Calculating solvent accessible surface area for single and combined chains  
"""

ids=[]
sasa=[]
id1=[]
id2=[]
comb=[]
expectedsum=[]
probability=[]
pro_sum=[]
totalSurfaceAreasComplex=[]
max_values=[]
chain_a=[]
chain_b=[]
max_dist={}
int_pro=[]
dist=[]

# for pdbs in pdb_id:
#     chains=glob.glob(f"pdb_new_files/{pdbs}/{pdbs}-*.pdb")        
#     for i in range(len(chains)):
#         structure = freesasa.Structure(chains[i]) #load the structure with index i   
#         result = freesasa.calc(structure) #calculate sasa for selected structure 
#         total=result.totalArea()
#         sasa.append(total)
#         ids.append(chains[i].split("/")[2].split(".")[0])

# # Creating dataframe for the SASA calculation
# # prep=list(zip(ids,sasa))
# # df=pd.DataFrame(prep,columns=["PDB ID","SASA"])


# print("single chain SASA calculation: √")

max_ppi_probs = []
max_ppi_pairs = []
i=0
#Estimating SASA for combined chains 
os.chdir("/Users/moustafashehata/Desktop/ppi/")
for pdbs in pdb_id:
    chains=glob.glob(f"pdb_new_files/{pdbs}/{pdbs}-*.pdb")       
    print(i, pdbs, chains)
    i+= 1
    
    max_prob = 0.0
    max_ppi_pair = (-1,-1)
    for i in range(len(chains)):
        for j in range( i+1, len(chains)):
            #print(chains[i],chains[j])
            
            os.system("cat %s %s > temp.pdb" % (chains[i], chains[j])) #combine two proteins and name it temp
            
            complexStructure = freesasa.Structure("temp.pdb")  # load complex structure 
            
            resultcomplex = freesasa.calc(complexStructure)    # calculate sasa for complex 
            
            ind1 = round((freesasa.calc(freesasa.Structure(f"{chains[i]}"))).totalArea(),4) #chain 1
            
            ind2 = round((freesasa.calc(freesasa.Structure(f"{chains[j]}"))).totalArea(),4) #chain 2
                        
            
            Sums = resultcomplex.totalArea()  #SUM
            
            expectedSum = ind1 + ind2  #expected sum 
            
            pro_sum.append(expectedSum)
                        
            probability =(expectedSum-Sums)/expectedSum
            
            int_pro.append((expectedSum-Sums)/expectedSum) #Interaction probability 
            
            chain_a.append(ind1) #SASA of chain 1 
            
            chain_b.append(ind2) #SASA of chain 2
            
            totalSurfaceAreasComplex.append(Sums)
            
            dist.append(probability)
                        
            id1.append(chains[i])
            
            id2.append(chains[j])
            
            comb.append(Sums)
            
#             if probability >= max_prob:
#                 max_prob = probability
#                 max_ppi_pair = (chains[i].split("/")[2].split(".")[0],chains[j].split("/")[2].split(".")[0])
                
#     max_ppi_probs.append(max_prob)
#     max_ppi_pairs.append([max_prob,pdbs,max_ppi_pair])
    
prep_comb=list(zip([text.split("/")[2].split(".")[0] for text in id1],
                   chain_a,
                   [text.split("/")[2].split(".")[0] for text in id2],
                   chain_b,
                   [round(num,2) for num in pro_sum], 
                   [round(num,2) for num in comb],
                   [round(num,4) for num in int_pro]))

df_comb=pd.DataFrame(prep_comb,columns=["PDB ID 1", "SASA 1" ,"PDB ID 2", "SASA 2", "Expected SASA sum", "SASA Sum","Interaction probability"])

df_comb.to_csv("preselection.csv",header=True,index=None)

df_comb[(df_comb["Interaction probability"] != -0.0000) | (df_comb["Interaction probability"] != 0.0000)].to_csv("selected_chains_for_docking.csv",header=True,index=None)

print("Combined chains SASA calculation: √")




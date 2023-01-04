import pandas as pd
import os

start = 0
end = 10

df = pd.read_csv("selected_chains.csv",sep=",")

#creat directories for PDB chains 
#for items, rows in df.iterrows():
#    pdb = rows["PDB ID"]
#    chain_a = rows["PDB ID 1"]
#    chain_b = rows["PDB ID 2"]

#    os.mkdir(f"../cport/{chain_a}")
#    os.system(f"cp {pdb}/{chain_a}.pdb ../cport/{chain_a}")

#    os.mkdir(f"../cport/{chain_b}")
#    os.system(f" cp {pdb}/{chain_b}.pdb ../cport/{chain_b}")

os.chdir("../cport")

while end <= 470:
    with open(f"../cport/cport_{start}_{end}.sh","w") as f:
        f.write("#!/bin/bash\n")
        df = pd.read_csv("selected_chains.csv",sep=",")
        df = df.iloc[start:end]
        for items, rows in df.iterrows():
            pdb = rows["PDB ID"]
            chain_a = rows["PDB ID 1"]
            chain_b = rows["PDB ID 2"]
            chain_a_id = chain_a.split("_")[1]
            chain_b_id = chain_b.split("_")[1]
            f.write(f"cd {chain_a}/\n")
            f.write(f"cport {chain_a}_clean.pdb {chain_a_id} --pred fast\n")
            f.write(f"cd ../\n")
            f.write(f"cd {chain_b}/\n")
            f.write(f"cport {chain_b}_clean.pdb {chain_b_id} --pred fast\n")
            f.write("cd ../\n")
    os.system(f"sbatch -e error_{start}_{end}.txt -o log_{start}_{end}.txt -t1440 -c 20 --mem=64G cport_{start}_{end}.sh -y {start}_{end}")
    start = end
    end += 10

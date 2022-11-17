#!/usr/bin/env python 
import pandas as pd
import os 

#Extract structure with RMSD value <= 3, 5, and 10 for 0.3 threshold dockng
rmsd_3_it0=[]
rmsd_5_it0=[]
rmsd_10_it0=[]
rmsd_5_3_it0=[]
rmsd_10_5_it0=[]

rmsd_3_it1=[]
rmsd_5_it1=[]
rmsd_10_it1=[]
rmsd_5_3_it1=[]
rmsd_10_5_it1=[]

rmsd_3_it1_water=[]
rmsd_5_it1_water=[]
rmsd_10_it1_water=[]
rmsd_5_3_it1_water=[]
rmsd_10_5_it1_water=[]

path=str(input("Please enter the folder path:\n"))
         
def extract_pdbs_3(x, y: float, z: str):
    for items,rows in x.iterrows():
        if rows["rmsd_all"] <= y:
            globals()[f"rmsd_{z}"].append(rows["#struc"])
            
def extract_pdbs_5(x, y: float, z: str):
    for items,rows in x.iterrows():
        if rows["rmsd_all"] >= y:
            globals()[f"rmsd_{z}"].append(rows["#struc"])
            
def extract_pdbs_10(x, y: float, z: str):
    for items,rows in x.iterrows():
        if rows["rmsd_all"] >= y:
            globals()[f"rmsd_{z}"].append(rows["#struc"])

def extract_pdbs_5_3(x, y: float, z: str):
    for items,rows in x.iterrows():
        if y > rows["rmsd_all"] >= 3:
            globals()[f"rmsd_{z}"].append(rows["#struc"])
    #print(globals()[f"rmsd_{z}"])
        
            
def extract_pdbs_10_5(x, y: float, z: str):
    for items,rows in x.iterrows():
        if y > rows["rmsd_all"] >= 5:
            globals()[f"rmsd_{z}"].append(rows["#struc"])
            
for file in os.listdir(f"analysis/{path}"):
    df_it0=pd.read_csv(f"analysis/{path}/{file}/it0/structures_rmsd-sorted.dat",sep="\s+")
    df_it1=pd.read_csv(f"analysis/{path}/{file}/it1/structures_rmsd-sorted.dat",sep="\s+")
    df_it1_water=pd.read_csv(f"analysis/{path}/{file}/it1/water/structures_rmsd-sorted.dat",sep="\s+")
    
    #Extract it0
    extract_pdbs_3(df_it0, 3 , "3_it0" )
    extract_pdbs_5(df_it0, 5 , "5_it0" )
    extract_pdbs_10(df_it0, 10 ,"10_it0")
    extract_pdbs_5_3(df_it0,5, "5_3_it0")
    extract_pdbs_10_5(df_it0,10,"10_5_it0")
    
    #Extract it1
    extract_pdbs_3(df_it1, 3, "3_it1")
    extract_pdbs_5(df_it1, 5, "5_it1")
    extract_pdbs_10(df_it1, 10, "10_it1")
    extract_pdbs_5_3(df_it1,5, "5_3_it1")
    extract_pdbs_10_5(df_it1,10,"10_5_it1")
    
    #Extract it1 water
    extract_pdbs_3(df_it1_water, 3, "3_it1_water")
    extract_pdbs_5(df_it1_water, 3, "5_it1_water")
    extract_pdbs_10(df_it1_water, 3, "10_it1_water")
    extract_pdbs_5_3(df_it1_water,5, "5_3_it1_water")
    extract_pdbs_10_5(df_it1_water,5, "10_5_it1_water")
    
print("")
print("*****IT0*****")
print("it0: <=3:",len(rmsd_3_it0))
print("it0: >= 5:",len(rmsd_5_it0))
print("it0: >= 10:",len(rmsd_10_it0))
print("it0: <5 x >=3:",len(rmsd_5_3_it0))
print("it0: <10 x >=5:",len(rmsd_10_5_it0))
print("*"*50)
print("*****IT1*****")
print("it1: <=3:",len(rmsd_3_it1))
print("it1: >= 5:",len(rmsd_5_it1))
print("it1: >= 10:",len(rmsd_10_it1))
print("it1: <5 x >=3:",len(rmsd_5_3_it1))
print("it1: <10 x >=5:",len(rmsd_10_5_it1))
print("*"*50)
print("*****IT1 water*****")
print("it1: <=3:",len(rmsd_3_it1_water))
print("it1: >= 5:",len(rmsd_5_it1_water))
print("it1: >= 10:",len(rmsd_10_it1_water))
print("it1: <5 x >=3:",len(rmsd_5_3_it1_water))
print("it1: <10 x >=5:",len(rmsd_10_5_it1_water))
print("*"*50)

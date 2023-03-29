import os
import csv
import glob
import pbxplore as pbx

def convert_pdb_pbx( pdb_file ):

    # Load the structure pdb files
    structure_reader = pbx.chains_from_files([pdb_file])
    chain_name, chain = next(structure_reader)
    dihedrals = chain.get_phi_psi_angles()
    pb_seq= pbx.assign(dihedrals)


    # Directory of the output
    output_dir = "/home/oem/Desktop/ALEX/Project/data/pbx_output/"

    [dirname,filename] = os.path.split( pdb_file )
    output_file = os.path.join(output_dir, filename+".csv")

    # Open and read the file
    with open(output_file, "w", newline='') as f:
        writer = csv.writer(f, delimiter= ",")
        writer.writerow([pb_seq])
        

    # Print confirmation message
    print("Resultados guardados en: {}".format(output_file))

if __name__ == "__main__":
    for pdb_path in glob.glob("/home/oem/Desktop/ALEX/Project/data/pdb_files/swissprot_pdb_v4/*.pdb"):
        convert_pdb_pbx( pdb_path )


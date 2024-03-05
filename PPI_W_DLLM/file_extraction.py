import tarfile
import os

def extract_pdb_names_from_tgz(tgz_file, output_file):
    with tarfile.open(tgz_file, "r:gz") as tar:
        pdb_files = [member for member in tar.getmembers() if os.path.splitext(member.name)[1] == '.pdb']
        with open(output_file, 'w') as f:
            for pdb_file in pdb_files:
                pdb_name = os.path.splitext(os.path.basename(pdb_file.name))[0]
                first_part, second_part = pdb_name.split('-')[:2]
                f.write(f"{first_part}\t{second_part}\n")

def process_tgz_files_in_directory(work_dir):
    tgz_files = [file for file in os.listdir(work_dir) if file.endswith(".tgz")]
    for tgz_file in tgz_files:
        output_file = os.path.splitext(tgz_file)[0] + ".txt"
        extract_pdb_names_from_tgz(os.path.join(work_dir, tgz_file), os.path.join(work_dir, output_file))
        print(f"Processed {tgz_file} and saved output to {output_file}")

#  workdir pathway !!!!
work_dir = "/home/pc550/Documents/PPI_W_DLLM/workdir"

# function to read all the .tgz files in a workdir
process_tgz_files_in_directory(work_dir)
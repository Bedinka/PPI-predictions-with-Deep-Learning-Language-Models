import gzip
import os
import shutil

def process_gz_files_in_directory(work_dir):
    for gz_file in os.listdir(work_dir):
        if gz_file.endswith(".gz"):
            gz_path = os.path.join(work_dir, gz_file)
            output_dir = '/home/dina/Documents/PPI_WDLLM/workdir/Alpha_pdb'
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, gz_file[:-3])  # Remove .gz extension
            with gzip.open(gz_path, 'rb') as f_in:
                with open(output_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print(f'Extracted file {gz_file} to {output_path}')

# workdir pathway
work_dir = "/home/dina/Documents/PPI_WDLLM/workdir/swiss/"

# function to read all the .gz files in a workdir
process_gz_files_in_directory(work_dir)

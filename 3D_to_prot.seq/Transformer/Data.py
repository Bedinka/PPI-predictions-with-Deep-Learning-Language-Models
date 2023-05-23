import os
import csv

csv_path = '/home/oem/Desktop/ALEX/Project/data/all_data.csv'
src_col_name = 'fasta_seq'
trg_col_name = '1d_seq'
output_dir = '/home/oem/Desktop/ALEX/Project/Transformer-master/data'
src_file_name = 'fasta.txt'
trg_file_name = 'structure.txt'

with open(csv_path, newline='') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)
    src_data, trg_data = zip(*[(row[2], row[1]) for row in reader])

with open(os.path.join(output_dir, src_file_name), 'w') as f:
    f.write('\n'.join(src_data))

with open(os.path.join(output_dir, trg_file_name), 'w') as f:
    f.write('\n'.join(trg_data))
import csv
import os 

def main(input_file, interactome_file, output_file, dir ):
    representative_sequences = set()
    input_path = os.path.join( dir, input_file )
    with open(interactome_file , 'r') as file:
        sequence_id = None
        for line in file:
            if line.startswith(">"):
                sequence_id = line.strip()[1:]
            else:
                sequence = line.strip()
                representative_sequences.add(sequence)
    print(f'Representative sequences:{len(representative_sequences)}')
    total_sequences = 0
    removed_sequences = 0
    interact_pos = 0
    interact_neg = 0

    with open(input_file, 'r', newline='', encoding='utf-8') as input_file:
        reader = csv.reader(input_file, delimiter='\t')
        header = next(reader)

        with open(output_file, 'w', newline='', encoding='utf-8') as output_file:
            writer = csv.writer(output_file, delimiter='\t')
            writer.writerow(header)

            for row in reader:
                print('Looking for Non rep in rows of input')
                total_sequences += 1
                if len(row) >= 3:  
                    prot_ID = row[0]
                    interact = row[1]
                    sequences = row[2]

                    seq1, seq2 = sequences.split('_')
                    seq1 = seq1.strip("[]'")
                    seq2 = seq2.strip("[]'")
                    representative_found = False
                    if seq1 and seq2 in representative_sequences:
                        representative_found = True
                        #print(f"Representative sequence was extracted {prot_ID}")
                    if representative_found:
                        writer.writerow(row)
                    else:
                        removed_sequences += 1
                        if interact == 0:
                            interact_neg += 1
                        else:
                            interact_pos +=1
                        print (f"Non-representative sequence was removed {prot_ID}")

    print("Redundancy removal: DONE ...")
    
    if total_sequences > 0:
        removal_percentage = (removed_sequences / total_sequences) * 100
        print(f"Total sequences processed: {total_sequences}")
        print(f"Sequences removed: {removed_sequences}")
        print(f"Percentage of sequences removed: {removal_percentage:.2f}%")
        print(f"Total NEGATIVE interact sequences removed: {interact_neg}")
        print(f"Total POSITIVE interact sequences removed: {interact_pos}")
    else:
        print("No sequences processed.")

if __name__ == "__main__":
    input_file = 'bert_train_9.tsv'
    interactome_file = 'interactom_nonredundant'
    output_file = 'bert_filtered_file.tsv'
    main(input_file, interactome_file, output_file, dir)

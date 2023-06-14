#!/usr/bin/env python3

from Bio import SeqIO
import shelve
import os
import csv

def read_seq( file ):
    sequences = SeqIO.parse(file, 'fasta')
    return {record.id: record for record in sequences}
    
def extract_seq(seq_dic, protname, l):
    r = ''
    s = seq_dic[protname].seq
    for i in l:
        r+=(s[i-1])
    return r

def extract_seqs(seq_dic, protname, ll):
    seq_list = []
    for l in ll:
        seq_list.append( extract_seq(seq_dic, protname, l) )
    return seq_list
           

def get_consecutive_motif(pos_list, threshold = 5):
    output = []
    consecutive_motif_list = []
    pos_list.sort()
    consecutive_motif = []
    prev_i = -3
    for pos in pos_list:
        dist = (pos - prev_i)
        if dist == 1:
            consecutive_motif.append( pos )
        elif dist == 2:
            consecutive_motif.append( pos-1 )
            consecutive_motif.append( pos )
        elif dist == 3:
            consecutive_motif.append( pos-2 )
            consecutive_motif.append( pos-1 )
            consecutive_motif.append( pos )
        else:
            if len(consecutive_motif) >= 1:
                consecutive_motif_list.append( consecutive_motif )
            consecutive_motif = [pos]
        prev_i = pos
    if len(consecutive_motif) >= 1:
        consecutive_motif_list.append( consecutive_motif )

    for consecutive_motif in consecutive_motif_list:
        if len(consecutive_motif) >= threshold:
            output.append( consecutive_motif )
    return output

def correctPN(proteinName):
    if proteinName[4] != ':':
        print('Protein name not changed:\t', proteinName)
        return proteinName
    print('Protein name changed:\t', proteinName[:4].upper()+proteinName[4:])
    return proteinName[:4].upper()+proteinName[4:]

class Data:
    def __init__(self, db):
        self.seq_dic = read_seq( 'multifasta.fa' )       
        self.db = shelve.open(db)
        print(len(self.db))
        self.domains = {}
        self.parse_domains()
        self.db.close()
    
    def parse_domains(self):
        for k in self.db:
            print(k)
            fields = k.split('__' )
            # if len(fields) != 4: continue
            [ proteinName, domainID, startPos, endPos ] = fields
            # set the protein name to the correct one, for the pdb database (XXXX:y), and for the alphafold database (AF-XXXXXX-F1-model_v4)
            proteinName = correctPN(proteinName)

            consecutive_motifs = get_consecutive_motif ( self.db[k] )
            
            # if len(consecutive_motifs) == 0: continue
            
            consecutive_motifs_seq = extract_seqs( self.seq_dic, proteinName, consecutive_motifs )
            
            info = [proteinName, int(startPos), int(endPos), self.db[k], consecutive_motifs, consecutive_motifs_seq ]
                        
            if domainID not in self.domains:
                self.domains[ domainID ] = [ info ]
            else:
                self.domains[ domainID ].append( info )


def find_species(proteinName):
    if proteinName[4] == ':':
        filepath = './data/pdb'+proteinName[:4].lower()+'.ent'
    else:
        filepath = './data/'+proteinName+'.pdb'
    good = ''
    with open(filepath, 'r') as f:
        for line in f:
            type_of_line = line.split()[0]
            if type_of_line == 'SOURCE':
                good += line
            elif type_of_line == 'REMARK':
                break
    good = good.replace('\t', '')
    good = good.replace('\n', '')
    good = good.replace(';', '')
    good = good.split()
    if 'ORGANISM_SCIENTIFIC:' in good:
        index = good.index('ORGANISM_SCIENTIFIC:')
        return good[index+1]+' '+good[index+2]
    elif 'ORGANISM_COMMON:' in good:
        index = good.index('ORGANISM_COMMON:')
        return good[index+1]+' '+good[index+2]
    elif 'SYNTHETIC:' in good:
        index = good.index('SYNTHETIC:')
        if good[index+1] == 'YES':
            return 'SYNTHETIC'
    else:
        return 'UNSEPCIFIED'

folder = os.getcwd().split('/')[-1]
print(folder)

domainInfo = shelve.open('domainInfo.shelve')

data = Data('domain_interactions.db')

table = open(folder+'_results.csv', 'w', encoding='utf-8')
writer = csv.writer(table)
writer.writerow(["Species", "PDB", "Domain", "Domain Name", "Motif seq", "Domain N-term", "Domain C-term", "Motif N-term", "Motif C-term", "Distance Domain-Motif"])

for domainID in data.domains:
    if len(data.domains[ domainID ]) == 1: continue
    print( "=======================" )
    print( domainID )
    domainName = domainInfo[domainID]
    for info in data.domains[ domainID ]:
        print( info )
        proteinName, startPos, endPos, interacting_residues, consecutive_motifs, consecutive_motifs_seq = info
        species = find_species(proteinName)
        for positions, seq in zip(consecutive_motifs, consecutive_motifs_seq):
            # distance = peptide position - domain position
            distance = min(positions)-endPos if max(positions)>startPos else max(positions)-startPos
            results = [species, proteinName, domainID, domainName, seq, startPos, endPos, min(positions), max(positions), distance]
            writer.writerow(results)

table.close()
print(folder)

# print( "=======================" )
# classify_peptides(['', 30, 80, [], [[22, 23], [400, 415, 420, 500]], ['THISSHOULDBECLOSE 22 23', 'THISSHOULDBEFAR 400 500']])
# print( "=======================" )
# classify_peptides_v2(['', 100, 180, [], [[22, 23], [70, 85], [200, 300], [400, 500]], ['Far from N-term 22, 23', 'Close to N-term 70, 85', 'Close to C-term 200, 300', 'Far from C-term 400, 500']])

def classify_peptides(info, threshold=20):
    proteinName, startPos, endPos, interacting_residues, consecutive_motifs, consecutive_motifs_seq = info
    index = 0
    for positions in consecutive_motifs:
        if max(positions)-startPos | max(positions)-startPos | min(positions)-endPos | min(positions)-endPos < threshold:
            print(consecutive_motifs_seq[index], 'is close to the domain')
        else:
            print(consecutive_motifs_seq[index], 'is far from the domain')
        index += 1

def classify_peptides_v2(info, threshold=20):
    proteinName, startPos, endPos, interacting_residues, consecutive_motifs, consecutive_motifs_seq = info
    for positions, seq in zip(consecutive_motifs, consecutive_motifs_seq):
        distance = min(positions)-endPos if max(positions)>startPos else max(positions)-startPos
        if abs(distance) < threshold:
            print(f'{seq} is close to the domains { "C-terminal" if max(positions)>startPos else "N-terminal"}')
        else:
            print(f'{seq} is far from the domains { "C-terminal" if max(positions)>startPos else "N-terminal"}')
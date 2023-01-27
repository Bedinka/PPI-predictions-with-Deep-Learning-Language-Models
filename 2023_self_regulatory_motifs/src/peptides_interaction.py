#!/usr/bin/env python3

from Bio import SeqIO
import shelve

def read_seq( file ):
    seq_dic = SeqIO.to_dict(SeqIO.parse(open(file), 'fasta'))
    return seq_dic
    
def extract_seq(seq_dic, protname, l):
    r = ''
    s = seq_dic[protname].seq
    for i in l:
        assert( i > 0 )
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
    prev_i = -1
    for pos in pos_list:
        if prev_i == pos - 1:
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
            #print(k)
            fields = k.split('__' )
            if len(fields) != 4: continue
            [ proteinName, domainID, startPos, endPos ] = fields  
            #print(self.db[k])
            consecutive_motifs = get_consecutive_motif ( self.db[k] )
                
            consecutive_motifs_seq = extract_seqs( self.seq_dic, proteinName, consecutive_motifs )
            #print(consecutive_motifs_seq)
            
            info = [proteinName, int(startPos), int(endPos), self.db[k], consecutive_motifs, consecutive_motifs_seq ]
            
            if len(consecutive_motifs) == 0: continue
            
            #print(domainID, info)
                        
            if domainID not in self.domains:
                self.domains[ domainID ] = [ info ]
            else:
                self.domains[ domainID ].append( info )

def classify_peptides(info, threshold=50):
    proteinName, startPos, endPos, interacting_residues, consecutive_motifs, consecutive_motifs_seq = info
    index = 0
    for positions in consecutive_motifs:
        if max(positions)-startPos | max(positions)-startPos | min(positions)-endPos | min(positions)-endPos > threshold:
            print(consecutive_motifs_seq[index], 'is close to the domain')
        else:
            print(consecutive_motifs_seq[index], 'is far from the domain')
        index += 1
            
            
data = Data('domain_interactions.db')


for domainID in data.domains:
    if len(data.domains[ domainID ]) == 1: continue
    print( "=======================" )
    print( domainID )
    for info in data.domains[ domainID ]:
        print( info )
        classify_peptides(info)

print( "=======================" )
classify_peptides(['', 30, 80, [], [[22, 23], [400, 415, 420, 500]], ['THISSHOULDBECLOSE', 'THISSHOULDBEFAR']])

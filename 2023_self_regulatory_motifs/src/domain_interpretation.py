#!/usr/bin/env python3

from Bio import SearchIO

def overlap(a, b):
    if a is None or b is None:
        return False
    return a[0] <= b[0] <= a[1] or b[0] <= a[0] <= b[1]

def domtbl_parser(file, E=10e-6):
    qrdict = {}
    evalue_filter = lambda hsp: hsp.evalue < E
    
    with open(file, 'r') as f:
        # transforms to zero-based and half-open intervals
        domtbl = SearchIO.parse(f, 'hmmscan3-domtab')
        for qresult in domtbl:
            print('\nPROTEIN:\t', qresult.id, '\n-------------------------\n')
            for hit in qresult.hits:
                # filter the hsp on the hit by evalue
                filtered_hit = hit.filter(evalue_filter)
                # if no hsps of the hit pass the threshold pop the hit from the QueryResult
                if filtered_hit is None:
                    qresult.pop(hit.id)
                else:
                    print(filtered_hit, '\n-----------\n', hit)
                    #save the "cleaned" hit to the QueryResult
                    qresult[hit.id] = filtered_hit
            # Only save QueryResults that contain information
            if len(qresult) > 0:
                print('ye \t\t\t\t', qresult.id)
                qrdict[qresult.id] = qresult
            else:
                print('no \t\t\t\t', qresult.id)
    return qrdict



def best_hits(qresult):
    for hit in qresult.hits:
        print(hit.id)
        print()
        for hsp in hit.hsps:
            print(hsp)
            print()
            print(hsp.evalue)
            print(hsp.query_range)


file = 'domtbl.out'

qrdict = domtbl_parser(file)

#best_hits(qrdict['AF-P49137-F1-model_v4'])
print('\n----------------------------------------\n')
# best_hits(qrdict['AF-A0A024RBG1-F1-model_v2'])

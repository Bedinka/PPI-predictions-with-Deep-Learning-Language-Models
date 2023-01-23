#!/usr/bin/env python

import sys
from Bio import SearchIO

def overlap(a,b):
    if a is None or b is None:
        return False
    return a[0] <= b[0] <= a[1] or b[0] <= a[0] <= b[1]

def domtbl_parser(file, E = 10e-6):
    qrdict = dict()
    evalue_filter = lambda hsp: hsp.evalue < E

    with open(file, 'r') as f:
        # transforms to zero-based and half-open intervals
        domtbl = SearchIO.parse(f, 'hmmscan3-domtab')
        for qresult in domtbl:
            for hit in qresult.hits:
                # filter the hsp on the hit by evalue
                filtered_hit = hit.filter(evalue_filter)
                # if no hsps of the hit pass the threshold pop the hit from the QueryResult
                if filtered_hit is None:
                    qresult.pop(hit.id)
                else:
                    prevhsp = None
                    index = 0
                    for hsp in filtered_hit:
                        # if two hsps overlap pop the worst one (they are sorted by evalue)
                        if overlap(prevhsp, hsp.query_range):
                            filtered_hit.pop(index)
                        prevhsp = hsp.query_range
                        index += 1
                    #save the "cleaned" hit to the QueryResult
                    qresult[hit.id] = filtered_hit
            # Only save QueryResults that contain information
            if len(qresult) > 0:
                qrdict[qresult.id] = qresult
        return qrdict


def best_hits(qresult):
    for h in qresult:
        for hsp in h:
            print(hsp)
            print()
            print(hsp.evalue)
            print(hsp.query_range)

file = sys.argv[1]

qrdict = domtbl_parser(file)

best_hits(qrdict['AF-A3KMH1-F1-model_v4'])
print('\n----------------------------------------\n')
best_hits(qrdict['AF-A0A024RBG1-F1-model_v2'])


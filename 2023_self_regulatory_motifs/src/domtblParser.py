#!/usr/bin/env python3

import shelve
from Bio import SearchIO

def overlap(a,b):
    if a is None or b is None:
        return False
    return a[0] <= b[0] <= a[1] or b[0] <= a[0] <= b[1]


def domtbl_parser(inputfile, outputfile, E = 10e-6):
    qrdict = shelve.open(outputfile)
    evalue_filter = lambda hsp: hsp.evalue < E

    with open(inputfile, 'r') as f:
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
                print(qresult.id)
                qrdict[qresult.id] = qresult
    qrdict.close()

domtbl_parser('domtbl.out', 'domdict.shelve')


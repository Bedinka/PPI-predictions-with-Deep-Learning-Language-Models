#!/usr/bin/env python3

from Bio import SwissProt
import csv

table = open('autoinhibitory_sprot.csv', 'w')
writer = csv.writer(table)
writer.writerow(["Accession", "Organism", "Type", "Location", "Sequence", "Note"])

for record in SwissProt.parse(open('uniprot_sprot.dat')):
    for i in record.features:
        if i.type == 'MUTAGEN': continue
        if 'note' in i.qualifiers.keys():
            if i.qualifiers['note'].lower().find('autoinhibitory') != -1:
                writer.writerow([
                    record.accessions[0],
                    record.organism,
                    i.type,
                    i.location,
                    i.location.extract(record.sequence),
                    i.qualifiers['note']
                ])
                print(record.accessions[0])
                print(record.organism)
                print(i.type)
                print(i.location)
                print(i.qualifiers['note'])
                print(i.location.extract(record.sequence))
                print('--------------\n')
table.close()

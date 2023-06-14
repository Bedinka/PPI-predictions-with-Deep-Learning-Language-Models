#!/usr/bin/env python3

import csv

def extact_uniprot_from_AF(af_id):
    s = af_id.split('-')
    if s[0] == 'AF':
        return s[1]
    return af_id

def overlap(a, b):
    if a is None or b is None:
        return False
    return a[0] <= b[0] <= a[1] or b[0] <= a[0] <= b[1]

table1 = open('SwissProt/SwissProt_results.csv', 'r')
reader1 = csv.reader(table1)
next(reader1)

table2 = open('autoinhibitory_sprot.csv', 'r')
reader2 = csv.reader(table2)
next(reader2)

table3 = open('sprot_found.csv', 'w')
writer1 = csv.writer(table3)
writer1.writerow(["Species","Accession","Domain","Domain Name","Motif seq","Domain N-term","Domain C-term","Motif N-term","Motif C-term","Distance Domain-Motif", "Curated motif Position", "Curated Seq"])

table4 = open('sprot_not_found.csv', 'w')
writer2 = csv.writer(table4)
writer2.writerow(["Species","Accession","Domain","Domain Name","Motif seq","Domain N-term","Domain C-term","Motif N-term","Motif C-term","Distance Domain-Motif", "Curated motif Position", "Curated Seq"])


motif_dic = {}
for line in reader1:
    id = extact_uniprot_from_AF(line[1])
    if id not in motif_dic.keys():
        motif_dic[id] = [line]
    else:
        motif_dic[id].append(line)

found = {}
# not_found = {}

for line in reader2:
    if line[0] in motif_dic.keys():
        for data in motif_dic[line[0]]:
            # if the positions overlap we will condsider that we found the motif/domain
            if overlap(tuple((data[7], data[8])), line[3].strip('[]').split(':')):
                writer1.writerow(data+[line[3], line[4]])
                print('\n', line[0], ':')
                print(data)
                print('-------------')
                print(line)
                if line[0] not in found.keys():
                    found[line[0]] = [data]
                else:
                    found[line[0]].append(data)
            # else:
            #     writer2.writerow(data+[line[3], line[4]])
            #     if line[0] not in not_found.keys():
            #         not_found[line[0]] = [data]
            #     else:
            #         not_found[line[0]].append(data)

print(len(found), 'found out of 289 (', 289-len(found), 'not found)')

table1.close()
table2.close()
table3.close()
table4.close()

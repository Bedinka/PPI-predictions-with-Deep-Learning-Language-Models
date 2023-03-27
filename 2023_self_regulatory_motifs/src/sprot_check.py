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

sprot_dic = {}
for line in reader2:
    if line[0] not in sprot_dic.keys():
        sprot_dic[line[0]] = [line[1:]]
    else:
        sprot_dic[line[0]].append(line[1:])

missing_sprot = {}

for line in reader1:
    id = extact_uniprot_from_AF(line[1])
    if id in sprot_dic.keys():
        new_list = []
        for data in sprot_dic[id]:
            # if the positions overlap we will condsider that we found the motif/domain
            if overlap(tuple((line[6], line[7])), data[2].strip('[]').split(':')):
                continue
            # if the motif found is a substring of the swissprot sequence we consider that we faound the motif/domain
            elif line[3] in data[3]:
                continue
            new_list.append(data)
        if new_list == []:
            sprot_dic.pop(id)
        else:
            sprot_dic[id] = new_list

c=0
for key,values in sprot_dic.items():
    print(key)
    c+=1
    for value in values:
        for i in value:
            print('\t', i)
        print()
    print('-----------------')
print(c)

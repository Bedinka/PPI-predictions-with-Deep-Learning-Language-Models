## main function ##
import csv
import os

os.path.abspath(os.path.curdir)
file_path = os.path.abspath(os.path.join(os.path.curdir, 'literature_pdb_residues_roles.csv'))
of = open(file_path, 'r')
reader = csv.reader(of)
next(reader, None)

temp1 = list()
temp2 = list()
map1 = dict()
map2 = dict()

for item in reader :
    temp1.append(item[2])
    temp2.append(item[9])

_temp1 = list(set(temp1))
_temp2 = list(set(temp2))

for i in range(len(_temp1)) :
    map1[_temp1[i]] = i

for i in range(len(_temp2)) :
    map2[_temp2[i]] = i

data = [[0 for j in range(len(_temp2))] for i in range(len(_temp1))]

of.seek(0)
next(reader, None)

for item in reader :
    data[map1[item[2]]][map2[item[9]]] += 1

for i in range(len(_temp2)) :
    for j in range(len(_temp1)) :
        result = "ROLE : {:50} RESIDUE : {:10} COUNT : {}".format(_temp2[i], _temp1[j], data[j][i])
        print(result)

#        
f = open("out.csv","w")
f.write('Role'+','+'Residue'+','+'Count'+'\n')
for i in range(len(_temp2)) :
    for j in range(len(_temp1)) :
        result = "ROLE : {:50} RESIDUE : {:10} COUNT : {}".format(_temp2[i], _temp1[j], data[j][i])
        f.write(format(_temp2[i])+','+format(_temp1[j])+','+format(data[j][i])+'\n')

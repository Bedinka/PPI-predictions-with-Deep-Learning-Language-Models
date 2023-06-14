#!/usr/bin/env python3

import os
import shelve

x = shelve.open('domdict.shelve')
y = shelve.open('domainInfo.shelve')
for i in x:
    for j in x[i]:
        print(j.accession)
        print(j.description)
        y[j.accession] = j.description

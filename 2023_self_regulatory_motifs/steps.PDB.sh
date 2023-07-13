#!/bin/bash

#cd $1
pwd
#../src/pdb2fasta.py > pdb2fasta.txt
echo 'finished pdb 2 fasta'
#hmmscan --cpu 2 --domtblout domtbl.out ../db/Pfam-A.hmm multifasta.fa > hmmscan.out
echo 'finished hmmscan'
../src/domtblParser.py > domtblParser.txt
echo 'domtbl parsed'
../src/AlphaFoldParser.PDB.py > AlphaFoldParser.txt
echo 'finished AlphaFoldParser'
../src/domainInfo.py > domainInfo.txt
echo 'finished domainInfo' $1
../src/peptides_interaction.py > peptides_interaction.txt
echo 'finished All'
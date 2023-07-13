#!/bin/bash

cd $1
pwd
../src/pdb2fasta.py > pdb2_fasta.txt
echo 'finished pdb 2 fasta' $1
hmmscan --domtblout domtbl.out ../db/Pfam-A.hmm multifasta.fa > hmmscan.out
echo 'finished hmmscan' $1
../src/domtblParser.py > domtblParser.txt
echo 'domtbl parsed' $1
../src/AlphaFoldParser.py > AlphaFoldParser.txt
echo 'finished AlphaFoldParser' $1
../src/domainInfo.py > domainInfo.txt
echo 'finished domainInfo' $1
../src/peptides_interaction.py > peptides_interaction.txt
echo 'finished All' $1

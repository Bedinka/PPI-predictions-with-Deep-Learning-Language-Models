#!/bin/bash

FILES="./fastas/*"

for f in $FILES
do
    hmmscan --domtblout domtbl.out ./db/Pfam-A.hmm $f > hmmscan.out &
    wait
    ./src/domain_interpretation.py domtbl.out
done
#!/bin/bash

i=0
for infile in fasta_clusters/*
do
    i=$((i + 1))
    echo $i / 351
    name=`echo $infile | cut -d . -f 1 | cut -d / -f 2`
    python clustalo.py --email bagren@cent.uw.edu.pl --stype protein --sequence $infile --outfile align/$name
done

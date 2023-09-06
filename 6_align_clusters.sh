#!/bin/bash


name=_representatives
if [ ! -f "align/$name.out.txt" ]
then
    python clustalo.py --email bagren@cent.uw.edu.pl --stype protein --sequence representatives_fasta.txt --outfile align/$name
fi
i=0
for infile in fasta_clusters/*
do
    i=$((i + 1))
    echo $i / 351
    name=`echo $infile | cut -d . -f 1 | cut -d / -f 2`
    if [ ! -f "align/$name.out.txt" ]
    then
        python clustalo.py --email bagren@cent.uw.edu.pl --stype protein --sequence $infile --outfile align/$name
    fi
done

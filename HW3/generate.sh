#!/bin/sh

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 LEN1 LEN2"
    exit 1
fi

LEN1=$1
LEN2=$2

python3 generate_fasta.py -seq1 $LEN1 -seq2 $LEN2
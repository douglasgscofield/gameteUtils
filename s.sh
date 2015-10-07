#!/bin/bash


for F in chr1.txt ; do
    echo "working on $F..."
    cut -f1,2 $F | /usr/bin/sed "s/^/$F:/" | tr ':' '\t' | head
done

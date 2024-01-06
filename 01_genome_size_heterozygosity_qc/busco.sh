#!/bin/bash

# genome assembly files: mazia_assembly.fna bedadeti_assembly.fna

for i in mazia bedadeti 
do
busco  \
-i ${i}_assembly.fna \
-m genome \
-l embryophyta_odb10 \
-c 16 \
-o busco.mazia.${i}.out \
--long \
-f
done 


fafafaf

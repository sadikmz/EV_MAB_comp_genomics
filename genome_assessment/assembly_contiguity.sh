#!/bin/bash

for i in mazia bedadeti
do 
~/apps/gfastats/build/bin/gfastats \
-b s \
-f EV_${i}.assembly.fna  \
-s s \
-t \
--seq-report \
--stats \
--sort descending > ${i}.gfastats.summmary.txt

done 

~/apps/gfastats/build/bin/gfastats \
-f purged.fa  \
-t \
--stats \
--cmd > summary_assembly_stat/assembly_${i}.hap1.stat
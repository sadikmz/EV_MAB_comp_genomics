#!/bin/bash

cpus=$cpus
genotype= #mazia or Bedadeit

#Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6754622/#S22

#1. Download and install Infernal 

# Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6754622/#S22
# cmpress Rfam.cm
#4.  Use cmscan to identify all subsequences that match any Rfam family with a score above the gathering cutoff (GA) selected by the Rfam curators:

~/apps/infernal-1.1.4/bin/cmscan --nohmmonly --rfam --cut_ga --fmt 2 --oclan --mpi --cpu 48 --oskip --clanin Rfam.clanin -o ${genotype}.cmscan.out --tblout ${genotype}.cmscan.tblout Rfam.cm ~/data/${genotype}.fna

#6. Fetch the high-scoring hits from the target sequence file 

esl-sfetch --index ~/data/${genotype}.fna

# 9. Using the .tblout file created in step 3, use esl-sfetch to fetch all subsequences that were a hit to an Rfam model, again using ‘awk’ and ‘grep’ as in step 6.

grep -v \^# ${genotype}.cmscan.tblout | awk ‘{ printf (“%s.%s/%d-%d %d %d %s\n”, $2, $4, $8, $9, $8, $9, $4); }’ | esl-sfetch -Cf ~/data/${genotype}.fna - > ${genotype}.rfam.fa
# The ${genotype}.rfam.fa file can now be used to map reads against to identify non-coding RNAs from families in Rfam. 
# Do not use Infernal (cmsearch or cmscan) to map the reads directly, as Infernal is not designed for short reads and will not give good results. 
# Use a dedicated read mapping program instead. Note that the fields (e.g. $2, $4) used in this command differ from those used in the similar command in step 6. 
# That is because the order and meaning of fields in the cmscan and cmsearch .tblout files differ.

# An alternative method for identifying ncRNAs from Rfam families in sequencing read datasets is to first map your reads against the complete target genome or dataset, and then look for overlaps between those mapped reads and Rfam hits found using cmscan (or cmsearch).

# 10. To do this, first convert the ${genotype}.cmscan .tblout file created in step 3 to GFF format:

grep -v ^\# ${genotype}.cmscan.tblout | awk ‘ { printf(“%s\tinfernal\t%s\t%s\t%s\t%s\t%s\t.\n”, $4, $2, $10, $11, $17, $12); }’ > ${genotype}.cmscan.gff
# The ${genotype}.cmscan.gff file can then be used as input to a program for identifying overlaps between two datasets, such as the bedtools intersect program (Quinlan, 2014).

# A similar command could be used to convert the cmsearch .tblout output from step 5 to GFF format
grep -v ^\# ${genotype}.cmsearch.tblout | awk ‘ { printf (“%s\tinfernal\t%s\t%s\t%s\t%s\t%s\t.\n”, $1, $3, $8, $9, $15, $10); }’ > ${genotype}.cmsearch.gff





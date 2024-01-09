#!/bin/bash


for i in mazia bedadeti
do 
# infernal 

cat infernal/${i}.cmscan.tblout | grep 'tRNA' > infernal/${i}.cmscan.tRNA.tblout 
# cat infernal/${i}.cmscan.tRNA.tblout | awk '{print $4,$11,$10}' OFS='\t' > infernal/${i}.cmscan.tRNA.bed 

awk '{if ($11 > $10) print $4"\t"$10"\t"$11; else print $4"\t"$11"\t"$10}' infernal/${i}.cmscan.tRNA.tblout | awk '{print $1"\t"$2"\t"$3}' OFS='\t'  > infernal/${i}.cmscan.tRNA.bed  

bedtools sort -i infernal/${i}.cmscan.tRNA.bed > infernal/${i}.cmscan.tRNA.sorted.bed

cat infernal/${i}.cmscan.tblout | grep -v '^#\|tRNA' > infernal/${i}.cmscan.non_tRNA.tblout 
done


# rRNA  - is bed file 

# non-coding RNA not capture in tRNAscan-SE but annotated in infernal  

bedtools intersect -b tRNAscan/out_mazia.bed -a infernal/mazia.cmscan.tRNA.sorted.bed -v -sorted > infernal/mazia.cmscan.non_overlapping_tRNA.bed
bedtools intersect -b tRNAscan/out_bedadeti.bed -a infernal/bedadeti.cmscan.tRNA.sorted.bed -v -sorted > infernal/bedadeti.cmscan.non_overlapping_tRNA.bed

# non-coding RNA not capture in tRNAscan-SE but homology search against A. thaliana and Rice in rRNA_homology_search.sh 

bedtools intersect -b tRNAscan/out_mazia.bed -a homology_search/mazia.rRNA.e05.blastn.out.sorted.merged.bed -v -sorted > homology_search/mazia.rRNA.e05.blastn.non_overlapping.bed
bedtools intersect -b tRNAscan/out_bedadeti.bed -a homology_search/bedadeti.rRNA.e05.blastn.out.sorted.merged.bed -v -sorted > homology_search/bedadeti.rRNA.e05.blastn.non_overlapping.bed

# non-coding RNA not capture in infernal but homology search against A. thaliana and Rice in rRNA_homology_search.sh 

bedtools intersect -b infernal/mazia.cmscan.non_overlapping_tRNA.bed -a homology_search/mazia.rRNA.e05.blastn.non_overlapping.bed -v -sorted > homology_search/mazia.rRNA.e05.blastn.non_overlapping.tRNAscan_infernal.bed
bedtools intersect -b infernal/bedadeti.cmscan.non_overlapping_tRNA.bed -a homology_search/bedadeti.rRNA.e05.blastn.non_overlapping.bed-v -sorted > homology_search/bedadeti.rRNA.e05.blastn.non_overlapping.tRNAscan_infernal.bed


# total tRNA

cat tRNAscan/out_mazia.bed infernal/mazia.cmscan.non_overlapping_tRNA.bed homology_search/mazia.rRNA.e05.blastn.non_overlapping.tRNAscan_infernal.bed > tRNAscan/mazia_infernal_tRNAscan_tRNA.bed 
cat tRNAscan/out_bedadeti.bed infernal/bedadeti.cmscan.non_overlapping_tRNA.bed homology_search/bedadeti.rRNA.e05.blastn.non_overlapping.tRNAscan_infernal.bed > tRNAscan/bedadeti_infernal_tRNAscan_tRNA.bed 



#!/bin/bash


cpus=$cpus
softmasked_genome=softmasked.fna # 
BASENAME=$(echo $softmasked_genome | sed "s/.fna//g")


# Prediction of protein coding genes using Braker

#1 BRAKER2 (v2.1.6) 

## installtaion: Use https://github.com/Gaius-Augustus/BRAKER or conda. 

conda create -n braker2 -c bioconda braker2

## Step 1. Generating hints from homologous proteins using ProtHint
## ProtHint installation: https://github.com/gatech-genemark/ProtHint 
## Proteins coding genes of 13 speceis downloaded from phytozome.v2: Arabidopsis thaliana, Oryza sativa, Oryza indica, Sorghum bicolor, Diascorea rotundata, 
## Brachypodium distachyon, Manihot esculenta, Zea mays, Ananas comosus, Musa acuminata, Musa balbisiana, Musa schizocarpa, and Ensete glaucum 

~/apps/ProtHint/bin/prothint.py your.repeatmasked_genome.fasta combined_multi_species_proteins.fa --threads $cpus

## activate braker's installation conda environment 
conda activate braker2  

## Run BRAKER2
braker.pl \
--species=ensete_ventricosum \
--genome=data/${BASENAME}.fna \
--hints=data/prothint_augustus.gff \
--bam=/data/${BASENAME}.rnaseq.bam \
--etpmode \
--softmasking \
--AUGUSTUS_ab_initio \
--cores $cpus 
#!/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --partition=hmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL,END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk

module load intel impi imkl

export OMP_NUM_THREADS=48

RNAseq=/home/lifesci/lfrwtp/data/rna_seq/DP8400009186BL_L01


# Prediction of protein coding genes using Braker

## RNA-seq alignment using hisat2 (v2.1.0)



hisat2-build  EV_mazia_v1.softmasked.1.fna EV_mazia_v1.softmasked.1.fna.index -p 48 --large-index

hisat2 \
-q -x EV_mazia_v1.softmasked.1.fna.index \
-p 48 \
-S EV_mazia.softmasked.1..rnaseq.sam \
-1 ${RNAseq}_573_1.fq.gz,${RNAseq}_574_1.fq.gz,${RNAseq}_575_1.fq.gz \
-2 ${RNAseq}_573_2.fq.gz,${RNAseq}_574_2.fq.gz,${RNAseq}_575_2.fq.gz \
--dta --dta-cufflinks \
--max-intronlen 160000 --no-mixed --very-sensitive --no-discordant --summary-file summary_file

#samtools view -@ 48 -bS mazia.v1.rnaseq.sam | samtools sort -@ 48 -o mazia.v1.rnaseq.sorted.bam 

samtools view -@ 48 -bS EV_mazia.softmasked.1.rnaseq.sam -o EV_mazia.softmasked.1.rnaseq.bam

samtools sort -@ 48 EV_mazia.softmasked.1.rnaseq.bam -o mazia.v1.rnaseq.sorted.bam

stringtie -o mazia_rnaseq.1.gff EV_mazia.softmasked.1.rnaseq.bam -p 48 -A gene_abundance.out


#1 BRAKER2 (v2.1.6) 

## installtaion: Use https://github.com/Gaius-Augustus/BRAKER or conda. Though compilation via the github documentation is the best way I would suggest to install the conda version.

conda create -n braker2 -c bioconda braker2

## Step 1. Generating hints from homologous proteins using ProtHint
## ProtHint installation: https://github.com/gatech-genemark/ProtHint 
## Proteins coding genes of 13 speceis downloaded from phytozome.v2: Arabidopsis thaliana, Oryza sativa, Oryza indica, Sorghum bicolor, Diascorea rotundata, 
## Brachypodium distachyon, Manihot esculenta, Zea mays, Ananas comosus, Musa acuminata, Musa balbisiana, Musa schizocarpa, and Ensete glaucum 

~/apps/ProtHint/bin/prothint.py your.repeatmasked_genome.fasta combined_multi_species_proteins.fa --threads 24

## activate braker's installation conda environment 
conda activate braker2  

## Run BRAKER2
braker.pl \
--species=ensete_ventricosum \
--genome=data/mazia_v1.500.softmasked.1.fna \
--hints=data/prothint_augustus.gff \
--bam=/data/mazia.v1.rnaseq.sorted.bam \
--etpmode \
--softmasking \
--AUGUSTUS_ab_initio \
--cores 48 \
--workingdir mazia.braker2.noUTR.out \
--makehub \
--geneMarkGtf mazia.braker2.noUTR.out/genemark.gtf \
--crf \
--email Sadik.Muzemil@warwick.ac.uk \
--nocleanup \
--MAKEHUB_PATH=/home/lifesci/lfrwtp/apps/BRAKER/MakeHub/  \
--GENEMARK_PATH=/home/lifesci/lfrwtp/apps/BRAKER/gmes_linux_64 


## add utr 

wd=mazia.braker2.noUTR.out

braker.pl \
--genome=data/mazia_v1.500.softmasked.1.fna \
--bam=/data/mazia.v1.rnaseq.sorted.bam \
--addUTR=on \
--flanking_DNA=149 \
--softmasking \
--workingdir=${wd} \
--AUGUSTUS_hints_preds=${wd}/augustus.hints.gtf \
--cores=48 \
--skipAllTraining \
--species=ensete_ventricosum01


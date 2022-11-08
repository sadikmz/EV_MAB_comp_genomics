#!/bin/bash
#SBATCH --job-name=trembd
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=3700
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL,END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk

module purge

module load GCC/10.2.0 CUDA/11.1.1 OpenMPI/4.0.5 

#mkdir fun_annot_out

# uniprot 
makeblastdb in ~/data/uniprot/uniprot_sprot.fasta -dbtype prot

# TrEMBL
gunzip ~data/uniprot_trembl.fasta.gz 
makeblastdb -in ~data/uniprot_trembl.fasta -dbtype prot 

# cog 
makeblastdb -in ~/data/cog_db/cog-20.fa -dbtype prot

# path to EggNOG diamond db
export EGGNOG_DATA_DIR=~/data/EggNOG


# Download eggNOG database
~/apps/eggnog/eggnog-mapper-2.1.8/download_eggnog_data.py -P -H -d 33090 -f --data_dir ~/data/EggNOG -y -q


# Run 
for i in mazia bedadeti
do 
# NCBI non-redundant db 

blastp \
-db ~/NCBI_NR/BLASTDB/nr \
-query mazia_bedadet_TE_excluded_AED25/${i}.25.TE_excluded.fasta \
-out fun_annot_out/${i}.maker2uni.cog.blastp.txt \
-evalue 0.000001 \
-outfmt 6 \
-num_alignments 1 \
-seg yes \
-soft_masking true \
-lcase_masking \
-max_hsps 1 \
-num_threads 36

# UNIPROT
blastp \
-db ~/data/uniprot/uniprot_sprot.fasta \
-query mazia_bedadet_TE_excluded_AED25/${i}.25.TE_excluded.fasta \
-out fun_annot_out/${i}.maker2uni.blastp.txt \
-evalue 0.000001 \
-outfmt 6 \
-num_alignments 1 \
-seg yes \
-soft_masking true \
-lcase_masking \
-max_hsps 1 \
-num_threads 48

# TeEMBL	
blastp \
-db ~data/uniprot_trembl.fasta \
-query mazia_bedadet_TE_excluded_AED25/${i}.25.TE_excluded.fasta \
-out fun_annot_out/${i}.maker2uni.trembl.blastp.txt \
-evalue 0.000001 \
-outfmt 6 \
-num_alignments 1 \
-seg yes \
-soft_masking true \
-lcase_masking \
-max_hsps 1 \
-num_threads 48

# cog

blastp \
-db ~/data/cog_db/cog-20.fa \
-query mazia_bedadet_TE_excluded_AED25/${i}.25.TE_excluded.fasta \
-out fun_annot_out/${i}.maker2uni.cog.blastp.txt \
-evalue 0.000001 \
-outfmt 6 \
-num_alignments 1 \
-seg yes \
-soft_masking true \
-lcase_masking \
-max_hsps 1 \
-num_threads 48


# Eggnog

emapper.py \
--cpu 48 \
-m diamond \
-o eggNOG \
-i /home/lifesci/lfrwtp/data/proteins/musa_ba.fasta \
--itype proteins \
--sensmode ultra-sensitive \
--output_dir eggnog_mapper_musa_ba.out  \
--dmnd_db /home/lifesci/lfrwtp/data/EggNOG/eggnog_proteins.dmnd \
--evalue 1e-5 \
--override \
--report_orthologs 

done 


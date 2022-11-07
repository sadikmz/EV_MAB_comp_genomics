#!/bin/bash 
#SBATCH --job-name=col
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

module purge

module load GCC/8.3.0 OpenMPI/3.1.4

## Collect annotation of gene models in fasta and gff files 

mkdir gene_models predicted_fasta 

gff3_merge -d mazia.maker.output/mazia_master_datastore_index.log -o gene_models/mazia.maker.gff
fasta_merge -d mazia.maker.output/mazia_master_datastore_index.log -o predicted_fasta/mazia
gff3_merge -g gene_models/mazia.maker.gff -o gene_models/mazia.maker.gene_models.gff

## filter gene models of annotation with minimum contig size by AED value 

### Collect splited annotation

for i in 0.16 0.21 0.26 0.31 0.36 0.41 0.46 0.51 0.56 0.61 0.66 0.71 0.76 0.81 0.86 0.91 0.96 1
   do
#     quality_filter.pl -a ${i} gene_models/mazia.maker.gene_models.gff >  gene_models/maker.${i}.gff
      grep -cP '\tgene\t' gene_models/maker.${i}.gff >> gene_utr_count/genes.AED.txt
      grep -cP '\t(three)_prime_UTR\t' gene_models/maker.${i}.gff >> gene_utr_count/utr3.AED.txt
      grep -cP '\t(five)_prime_UTR\t' gene_models/maker.${i}.gff >> gene_utr_count/utr5.AED.txt

#      paste AED_list.txt gene_utr_count/genes.AED.txt gene_utr_count/3utr.AED.txt gene_utr_count/5utr.AED.txt > gene_utr_count/genes.35utr.AED.txt
   done 

paste AED_list.txt gene_utr_count/genes.AED.txt gene_utr_count/utr3.AED.txt gene_utr_count/utr5.AED.txt > gene_utr_count/genes.35utr.txt
# merge count of genes and UTR into single file


#  Generate frequency table of gene models vs AED value 

AED_cdf_generator.pl -b 0.01 gene_models/bedadeti.maker.gene_models.gff > AED_unfiltred.genemodels.txt 


## filter gene models by their AED value 
## script source: https://www.biostars.org/p/364069/
##seqtk seq -l0 mazia.maker.proteins.0-0.30_AED.fasta | grep -v -e "AED:0.[3][5]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.30_AED.f.fasta

mkdir AED_fasta/AED_fasta_split
pred_proteins_fasta=predicted_fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[0-1][0-9]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.19_AED.fasta
seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[2][0]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.20_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.19_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.20_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.20_AED.fasta


seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[2][0-1]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.2-0.21_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.19_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.2-0.21_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.21_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[2][0-2]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.2-0.22_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.19_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.2-0.22_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.22_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[2][0-3]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.2-0.23_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.19_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.2-0.23_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.23_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[2][0-4]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.2-0.24_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.19_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.2-0.24_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.24_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[2][1-5]" |tr '\t' '\n' > AED_fasta/mazia.maker.proteins.0.21-0.25_AED.fasta
cat AED_fasta/mazia.maker.proteins.0.21-0.25_AED.fasta  AED_fasta/mazia.maker.proteins.0-0.20_AED.fasta > AED_fasta/mazia.maker.proteins.0-0.25_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[2][0-6]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.2-0.26_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.19_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.2-0.26_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.26_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[2][0-7]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.2-0.27_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.19_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.2-0.27_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.27_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[2][0-8]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.2-0.28_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.19_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.2-0.28_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.28_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[2][0-9]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.2-0.29_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.19_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.2-0.29_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.29_AED.fasta


seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[0-2][0-9]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.29_AED.fasta
seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[3][0]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.30_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.29_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.30_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.30_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[3][0-1]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.3-0.31_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.29_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.3-0.32_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.32_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[3][0-2]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.3-0.32_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.29_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.3-0.32_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.32_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[3][0-3]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.3-0.33_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.29_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.3-0.33_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.33_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[3][0-4]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.3-0.34_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.29_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.3-0.34_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.34_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[3][1-5]" |tr '\t' '\n' > AED_fasta/mazia.maker.proteins.0.31-35_AED.fasta
cat AED_fasta/mazia.maker.proteins.0.31-35_AED.fasta  AED_fasta/mazia.maker.proteins.0-0.30_AED.fasta > AED_fasta/mazia.maker.proteins.0-0.35_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[3][0-6]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.3-0.36_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.29_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.3-0.36_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.36_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[3][0-7]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.3-0.37_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.29_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.3-0.37_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.37_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[3][0-8]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.3-0.38_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.29_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.3-0.38_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.38_AED.fasta

seqtk seq -l0 ${pred_proteins_fasta}/mazia.all.maker.proteins.fasta| paste - - | grep -e " AED:0.[3][0-9]" |tr '\t' '\n' > AED_fasta/AED_fasta_split/mazia.maker.proteins.0.3-0.39_AED.fasta
cat AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.29_AED.fasta AED_fasta/AED_fasta_split/mazia.maker.proteins.0.3-0.39_AED.fasta > AED_fasta/AED_fasta_split/mazia.maker.proteins.0-0.39_AED.fasta


# run busco 

AED_list=$(cat AED_list.0.txt) # max AED values of splitted gene models (eg. 0-0.39_AED in mazia.maker.proteins.0-0.39_AED.fasta) 

#mkdir collect_split
AED_list="0-0.30_AED"
for i in ${AED_list}
do  

busco  -i ../AED_fasta/AED_fasta_split/mazia.maker.proteins.${i}.fasta \
-m protein \
-l embryophyta_odb10 \
-c 48 \
-o busco.${i}.out \
--long \
-f \
--offline 

cp busco.${i}.out/short_summary.specific.embryophyta_odb10.busco.${i}_AED.out.txt  collect_split

done 

cd collect_split

generate_plot.py -wd .


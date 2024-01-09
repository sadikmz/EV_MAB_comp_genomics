#!/bin/bash

# This script was used to extract mapping coverage of reads across coding sequences in enset and Musa species genomes. The input indicated here were mapping pan-Ensete ventricosum NGS reads against Musa spp. genomes
# to extract coverge reads across CDS. The same scrit was used to calculated CDS coverage of A-genome and B-genome NGS reads against enset genomes.
module purge
module load GCC/10.2.0 OpenMPI/4.0.5 GCCcore/10.3.0


#1.  mapping panEV reads against Musa acuminata and Musa balbisiana genomes. Prefix of PanEV reads were save as EV_WGS.list


PATH_GENOMES=/home/l/lfrwtp/data/genomes/
PATH_WGS=/home/l/lfrwtp/data/EV_wgs
cpus=48

## index genomes 
bwa-mem2 index $PATH_GENOMES/Musa_acuminata.fna
bwa-mem2 index $PATH_GENOMES/Musa_balbisiana.fna


for READS in $(cat EV_WGS.list)
do
query_genotype=$READS
ref_genome=Musa_acuminata
Fread=${PATH_WGS}/${query_genotype}.1_val_1.fq.gz
Rread=${PATH_WGS}/${query_genotype}.2_val_2.fq.gz
ref=${PATH_GENOMES}/${ref_genome}.fna
output_path=panEV_mapping_alignment
mkdir $output_path
mkdir ${output_path}/panEV_${ref_genome} 
output_dir=${output_path}/panEV_${ref_genome}
out=${query_genotype}_${ref_genome}
prefix=$out

## map PE reads  
## Run bwa-mem2
bwa-mem2 mem -t $cpus $$ref $Fread $Rread |  samtools view - -Sb -@$cpus | samtools view -b -@$cpus -F 4 -o ${output_dir}/$out.allMapped.bam

query_genotype=$READS
ref_genome=Musa_balbisiana
Fread=${PATH_WGS}/${query_genotype}.1_val_1.fq.gz
Rread=${PATH_WGS}/${query_genotype}.2_val_2.fq.gz
ref=${PATH_GENOMES}/${ref_genome}.fna
output_path=panEV_mapping_alignment
mkdir $output_path
mkdir ${output_path}/panEV_${ref_genome} 
output_dir=${output_path}/panEV_${ref_genome}
out=${query_genotype}_${ref_genome}
prefix=$out


## map PE reads  
## Run bwa-mem2
bwa-mem2 mem -t $cpus $$ref $Fread $Rread |  samtools view - -Sb -@$cpus | samtools view -b -@$cpus -F 4 -o ${output_dir}/$out.allMapped.bam

done 

## merge bams 
samtools merge panEV_mapping_alignment/panEV_Musa_acuminata.all.bam panEV_mapping_alignment/panEV_Musa_acuminata/*_Musa_acuminata.bam -@$cpus

## Remove PCR duplicates 

java -Xmx50G -jar $PICARD/picard.jar MarkDuplicates \
INPUT=$panEV_mapping_alignment/panEV_Musa_acuminata.all.bam \
O=panEV_mapping_alignment/panEV_Musa_acuminata.all.markdup.bam  \
M=panEV_mapping_alignment/panEV_Musa_acuminata.all.markdup.bam.metrics.txt \
REMOVE_DUPLICATES=True

java -Xmx50G -jar $PICARD/picard.jar MarkDuplicates \
INPUT=panEV_mapping_alignment/panEV_Musa_balbisiana.all.bam \
O=panEV_mapping_alignment/panEV_Musa_balbisiana.all.markdup.bam  \
M=panEV_mapping_alignment/panEV_Musa_balbisiana.all.markdup.bam.metrics.txt \
REMOVE_DUPLICATES=True

samtools view -b -@$cpus panEV_mapping_alignment/panEV_Musa_acuminata.all.markdup.bam | samtools sort - -@$cpus  -o panEV_mapping_alignment/panEV_Musa_acuminata.allMapped.sorted.bam
rm panEV_mapping_alignment/panEV_Musa_acuminata/*_Musa_acuminata.bam

samtools merge panEV_mapping_alignment/panEV_Musa_balbisiana.all.markdup.bam panEV_mapping_alignment/panEV_musa_acuminata/*_Musa_balbisiana.bam -@$cpus
samtools view -b -@$cpus panEV_mapping_alignment/panEV_Musa_balbisiana.bam | samtools sort - -@$cpus  -o panEV_mapping_alignment/panEV_Musa_balbisiana.allMapped.sorted.bam
rm panEV_mapping_alignment/panEV_Musa_balbisiana/*_Musa_balbisiana.bam

#2. Extract reads spanning coding sequences 
ref_genome=Musa_acuminata # and Musa_balbisiana 
out=panEV_mapping_alignment/panEV_${ref_genome}

## generate bed files gene and their complements (non-coding regions)

gff2bed < $PATH_GENOMES/${ref_genome}_maker_annotation.gff3 | grep 'gene' | cut -f1-3 > $PATH_GENOMES/${ref_genome}.gene.bed
samtools faidx $PATH_GENOMES/${ref_genome}.fna
cut -f1,2 $PATH_GENOMES/${ref_genome}.fna.fai > $PATH_GENOMES/${ref_genome}.sizes
bedtools complement -i $PATH_GENOMES/${ref_genome}.gene.bed -g $PATH_GENOMES/${ref_genome}.sizes > $PATH_GENOMES/${ref_genome}.gene.complement.bed 
cat $PATH_GENOMES/${ref_genome}.gene.bed $PATH_GENOMES/${ref_genome}.gene.complement.bed | sort | uniq > $PATH_GENOMES/${ref_genome}.gene_complement.bed 

## For gene and gene-complement 
samtools view -b -h -L $PATH_GENOMES/${query_genotype}.gene_complement.bed $out.allMapped.sorted.markdup.bam > $out.allMapped.gene_complement.bam
# samtools view -b -h -L $PATH_GENOMES/EV_${ref_genome}.gene_complement.bed $out.merged.primary.allMapped.bam > $out.primary.allMapped.gene_complement.bam 

## For gene only 
samtools view -b -h -L $PATH_GENOMES/${query_genotype}.gene.bed $out.allMapped.gene_complement.bam > $out.allMapped.gene.bam
# samtools view -b -h -L $PATH_GENOMES/EV_${ref_genome}.gene.bed $out.primary.allMals -pped.gene_complement.bam > $out.allMapped.primary.gene.bam

## index bam files 
samtools index $out.allMapped.gene.bam

## Coverage 
## bam coverage to bigwib: https://www.biostars.org/p/150036/
~/apps/T2T-Polish/coverage/sam2paf.sh $out.allMapped.gene.bam $out.allMapped.gene.paf $out.allMapped.gene

bamCoverage -b $out.allMapped.gene.bam --numberOfProcessors $cpus -o $out.allMapped.gene.coverage.bw
bedtools coverage -b $out.allMapped.gene.bam -a $PATH_GENOMES/${query_genotype}.gene.bed > $out.allMapped.gene.cov.bed



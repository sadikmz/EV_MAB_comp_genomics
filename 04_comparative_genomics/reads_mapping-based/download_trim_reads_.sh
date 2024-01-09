#!/bin/bash


```bash
cpus=48

#download Ensete ventricosum reads

# EV_landrace SRR files in PRJNA344540 were save as EV_SRR.txt  

esearch -db sra -query PRJNA344540 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | grep -f EV_SRR.txt | xargs -n 1 -P $cpus fastq-dump --split-files --gzip --skip-technical
esearch -db sra -query PRJNA252658 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | grep -E 'SRR1515268|SRR1515269' | xargs -n 1 -P $cpus fastq-dump --split-files --gzip --skip-technical

# prefix of all downloaded fastq files including the newlly generated NGS reads of Mazia were saved as EV_sra_list.txt)
ls -lht *.fastq.gz | awk '{print $9' | sed 's/_?.fastq.gz//g' > EV_sra_list.txt

# Trim adapter sequences and remove poor quality reads 

for i in $(cat EV_sra_list.txt)
do 
trim_galore -cores $cpus --quality 30 --fastqc --output_dir trim_out  --paired ${i}_?.fastq.gz 
done 
 
# for A-genome and B-genome reads of Musa species 
esearch -db sra -query PRJNA540118 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | grep -E 'SRR8989633|SRR9850640|SRR9734079|SRR8989632|SRR8989638|SRR8989628|SRR8989629|SRR8989634|SRR8989630|SRR9734077|SRR9850639|SRR9734074|SRR9734076' | xargs -n 1 -P $cpus fastq-dump --split-files --gzip --skip-technical
esearch -db sra -query PRJEB33317 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | grep 'ERR3412984' | xargs -n 1 -P $cpus fastq-dump --split-files --gzip --skip-technical
esearch -db sra -query PRJEB35002 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | grep -E 'ERR3606950|ERR3606951' | xargs -n 1 -P $cpus fastq-dump --split-files --gzip --skip-technical
esearch -db sra -query PRJEB58004 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | grep -E 'ERR10695603|ERR10695605|ERR10695608|ERR10695529|ERR10695565' | xargs -n 1 -P $cpus fastq-dump --split-files --gzip --skip-technical
esearch -db sra -query PRJNA216985 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | grep -E 'SRR956987|SRR956987' | xargs -n 1 -P $cpus fastq-dump --split-files --gzip --skip-technical
esearch -db sra -query PRJNA413600 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | grep -E 'SRR956987|SRR957627|SRR6147592' | xargs -n 1 -P $cpus fastq-dump --split-files --gzip --skip-technical
esearch -db sra -query PRJNA432894 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | grep -E 'SRR6996493|SRR6996492|SRR6996490|SRR6996491|SRR6996489|SRR6996488' | xargs -n 1 -P $cpus fastq-dump --split-files --gzip --skip-technical
esearch -db sra -query PRJNA540118 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | grep -E 'SRR8989633|SRR9850640|SRR9734079' | xargs -n 1 -P $cpus fastq-dump --split-files --gzip --skip-technical
esearch -db sra -query PRJNA540118 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | grep -E 'SRR8989633|SRR9850640|SRR9734079' | xargs -n 1 -P $cpus fastq-dump --split-files --gzip --skip-technical

# prefix of all downloaded fastq files A- and B-subgenomes were saved as EV_sra_list.txt)
ls -lht *.fastq.gz | awk '{print $9' | sed 's/_?.fastq.gz//g' > musaSpp_sra_list.txt

# adapter and poor quality reads trimming 
for i in $(cat musaSpp_sra_list.txt)
do 
trim_galore -cores 48 --quality 30 --fastqc --output_dir trim_out --paired ${i}_?.fastq.gz 
done
```bash
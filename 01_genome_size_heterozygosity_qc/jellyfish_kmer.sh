#!/bin/bash 
prog=~/apps/jellyfish/

# Trim fasta files  
trim_galore --cores 40 --quality 30 --fastqc -out jungle_seed --paired mazia_1.fastq.gz mazia_2.fastq.gz
trim_galore --cores 40 --quality 30 --fastqc -out jungle_seed --paired bedadeti_1.fastq.gz bedadeti_2.fastq.gz

ls -lht *trimmed.fq | awk '{print $9}' | sed 's/_._trimmed.fq.gz//g' | uniq > accession_list

# Run jellyfish and genomescope2
for K in 17 21 27 31;
do
# This will generate kmer count and kmer histogram for esimated genome size of 580M and with estimated 43X coverage.
    for L in $(cat accession_list)
        do
        # generate Kmer tabel 
        ${prog}/jellyfish-linux count -C -m $k -s 580M --bf-size 25G -t 16 <(zcat ${K}.${L}_1_trimmed.fq.gz) <(zcat ${K}${L}_2_trimmed.fq.gz) -o ${L}.${K}_mer.reads.jf
        
        # generate Kmer histogram
        ${prog}/jellyfish-linux histo -t 16 ${L}.${K}_mer.reads.jf > ${L}.${K}_mer.histo
        rm -rf ${L}.${K}_mer.reads.jf

        # Run genomescope 
        genomescope2 \
        - ${L}${K}_mer.histo \
        -o genomescope_${i}_k17 \
        -p 2 \
        -k ${K} \
        --testing \
        --fitted_hist
        done
done


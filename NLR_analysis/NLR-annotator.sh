#!/bin/bash

# Path to files and programs
genomes=path-to-ref-genomes # musa_ac.fna, musa_ba.fna, mazia.fna, bedadet.fna 
# 
prog=path-to-/NLR-annotator

mkdir NLR-annotator.out

for i in mazia bedadeti musa_ac musa_ba
do 
        java -Xmx100G -jar ${prog}/ChopSequence.jar \
        -i ${i}.fna \
        -o NLR-annotator.out/${i}_choseq_out.fasta

        java -Xmx100G -jar ${prog}/NLR-Parser.jar \
        -t 48 \
        -y ${prog}/meme_4.9.0/src/mast \
        -x ${prog}/meme.xml \
        -i NLR-annotator.out/${i}_choseq_out.fasta \
        -c NLR-annotator.out/${i}.00.nlr.xml

        java -Xmx100G -jar ${prog}/NLR-Annotator.jar \
        -i NLR-annotator.out/${i}.00.nlr.xml \
        -o NLR-annotator.out/${i}.nlr.txt  \
        -g NLR-annotator.out/${i}.NLR-annotator.out.gff \
        -b NLR-annotator.out/${i}.NLR-annotator.out.bed \
        -m NLR-annotator.out/${i}.NLR-annotator.out.motifs.bed \
        -a NLR-annotator.out/${i}.NLR-anotator.out.nbarkMotifAlignment.fasta 

        # extract predicted NLR fasta from bed file 

        bedtools getfasta \
        -fi ${i}.fna -bed NLR-annotator.out/${i}.NLR-annotator.out.bed > NLR-annotator.out/${i}.NLR-annotator.fasta 

        # Annotate NLR class 
        ./annotateClasses.v1.py ${i}.nlr.txt > ${i}.NLR_annotator_NLR_class.txt

done 

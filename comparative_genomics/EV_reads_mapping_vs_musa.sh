#!/bin/bash
#SBATCH --job-name=mapping
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL,END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk


set -eu
 
PATH_WGS=path-to-WGS-reads
musa_genomes=path-to-ref-genomes
PICARD=path-to-picard-2.27.1-0

# sed -E "s/^((@|\+)SRR[^.]+\.[^.]+)\.(1|2)/\1/" SRR3239806.fastq > SRR3239806.fixed.fastq

mkdir align_out 

ref_genomes=$(cat ref_genomes_list.txt) # ref genomes: musa_ac.fna musa_ba.fna 

for i in musa_ac musa_ba
    do 
        for j in mazia bedadeti # mazia_1_val_1.fq.gz, mazia_1_val_2.fq.gz  bedadeti_1_val_1.fq.gz bedadeti_2_val_2.fq.gz 
        do 
                bwa-mem2 mem -M -t 48 ${musa_genomes}/${i}.fna.index $PATH_WGS/${j}_1_val_1.fq.gz $PATH_WGS/${j}_2_val_2.fq.gz -o align_out/${i}.${j}.sam

                # convert sort sam

                java -Xmx50G -jar $PICARD/picard.jar SortSam \
                I=align_out/${i}.${j}.sam \
                O=align_out/${i}.${j}.sorted.sam \
                SORT_ORDER=coordinate \
                CREATE_INDEX=true \
                TMP_DIR=./tmp/
                VALIDATION_STRINGENCY=SILENT

                # convert sorted sam to bam format

                java -Xmx50G -jar $PICARD/picard.jar SamFormatConverter \
                I=align_out/${i}.${j}.sorted.sam \
                O=align_out/${i}.${j}.sorted.bam

                # Indexting bam files 

                java -Xmx50G -jar $PICARD/picard.jar BuildBamIndex \
                I=align_out/${i}.${j}.sorted.bam

                #calcualte alignment summary 

                java -Xmx50G -jar $PICARD/picard.jar CollectAlignmentSummaryMetrics \
                R=${musa_genomes}/Musa_acuminata_pahang_v4.fna \
                I=align_out/${i}.${j}.sorted.bam \
                O=align_out/${i}.${j}.sorted.bam.AlignMet.output.txt

                # # REMOVE_DUPLICATES mazia

                java -Xmx50G -jar $PICARD/picard.jar MarkDuplicates \
                INPUT=align_out/${i}.${j}.sorted.dedup.bam \
                O=align_out/${i}.${j}.markdup.merged.sorted.bam \
                M=align_out/${i}.${j}.markdup.merged.sorted.metrics.txt \
                REMOVE_DUPLICATES=True

                # # extract unligned sequence 
                # samtools view -f 12 align_out/${i}.${j}.markdup.merged.sorted.bam -o align_out/${i}.${j}.unaligned.bam
                # samtools bam2fq -1 $PATH_WGS/panreads/${i}/${i}.${j}.unaligned.11.fq -2 $PATH_WGS/panreads/${i}/${i}.${j}.unaligned.22.fq -n align_out/${i}.${j}.unaligned.bam -@ 36

                # mapping statistics 
        #        qualimap bamqc -bam align_out/${i}.${j}.sorted.bam -outdir align_out -outfile ${i}.${j}.sorted.bam.qualimap -sd --java-mem-size=10G

                # coverage 
                bedtools genomecov -ibam align_out/${i}.${j}.sorted.dedup.bam -bga > align_out/${i}.${j}.coverage.bed 

                # Remove created sam and bam files 

                rm align_out/*.sam
                rm align_out/*.bam
        done 

        
done



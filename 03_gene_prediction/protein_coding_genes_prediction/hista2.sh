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
cpus=$cpus
softmasked_genome=softmasked.fna # 
BASENAME=$(echo $softmasked_genome | sed "s/.fna//g")

## RNA-seq alignment using hisat2 (v2.1.0)

hisat2-build  $softmasked_genome ${softmasked_genome}.index -p $cpus --large-index

hisat2 \
-q -x ${softmasked_genome}.index \
-p $cpus \
-1 ${RNAseq}_573_1.fq.gz,${RNAseq}_574_1.fq.gz,${RNAseq}_575_1.fq.gz \
-2 ${RNAseq}_573_2.fq.gz,${RNAseq}_574_2.fq.gz,${RNAseq}_575_2.fq.gz \
--dta --dta-cufflinks \
--max-intronlen 160000 \
--no-mixed \
--very-sensitive \
--no-discordant \
--summary-file summary_file | samtools view -@$cpus -bS  | samtools sort -@$cpus -o ${BASENAME}.rnaseq.bam -T hisat2_tmp


stringtie -o ${BASENAME}_rnaseq.gff ${BASENAME}.rnaseq.bam  -p $cpus -A gene_abundance.out


#!/bin/bash
#SBATCH --job-name=trna1
#SBATCH --partition=cnode
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem-per-cpu=4571
#SBATCH --time=48:00:00
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk

MY_NUM_THREADS=$SLURM_CPUS_PER_TASK

export OMP_NUM_THREADS=$MY_NUM_THREAD

cpus=$cpus
data_dir=~/data/genomes.fna # 
BASENAME=$(echo $softmasked_genome | sed "s/.fna//g")
genotype=
## tRNAscan-SE (v2.0.9)

tRNAscan-SE -E ~/data/genotype.fna -HQ -X -L -y --threads $cpus -o# -f# -m# -s# -a# -b#--detail -p



## Ribosomal RNA (rRNA) fragments were predicted by BLASTN (e-value ≤ 1e−5) homology search of Arabidopsis (5S, 5.8S. and 18S rRNA) and Rice (28S rRNA) template rRNA sequences against the assembled genome.


genome_list=$(cat epo_mazia_bedadeti.genomes_list.txt)

makeblastdb -in ~/data/rRNA/rRNA.fna -input_type fasta -dbtype nucl


for i in $(ls -lht ~/data/*.fna| awk '{print $9}')
do
	BASENAME=$(basename $i | sed 's/.fna//g')

	blastn \
	-db rRNA.fna \
	-query ~/data/${BASENAME}.fna \
	-outfmt 6 \
	-evalue 1e-5 \
	-num_threads 48 \
	-num_alignments 1 \
	-out ${BASENAME}.rRNA.e05.blastn.out 

	# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore


	# convert blastn to bed 
	# swap bed coordinates that start position > end
	awk '{if ($7 > $8) print $1"\t"$8"\t"$7"\t"$2; else print $1"\t"$7"\t"$8"\t"$2}' ${BASENAME}.rRNA.e05.blastn.out > ${BASENAME}.rRNA.e05.blastn.out.bed

	#sort bed files
	bedtools sort -i ${BASENAME}.rRNA.e05.blastn.out.bed > ${BASENAME}.rRNA.e05.blastn.out.sorted.bed
	  
	# formatting bed file in case the new line character the end of each line  
	cat ${BASENAME}.rRNA.e05.blastn.out.sorted.bed | awk '{print $1, $2, $3, $4}' OFS='\t' > ${BASENAME}.rRNA.e05.blastn.out.sorted.1.bed
	bedtools merge -i ${BASENAME}.rRNA.e05.blastn.out.sorted.1.bed -c 4 -o distinct > ${BASENAME}.rRNA.e05.blastn.out.sorted.merged.bed

done 



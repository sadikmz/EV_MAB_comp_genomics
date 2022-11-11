#!/bin/bash
#SBATCH --job-name=NLR
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL,END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk

module purge

module load intel impi imkl

module load GCC/9.3.0
module load OpenMPI/4.0.3

#dmtcp_launch -i 3600 ./hmmsearch_script

#dmtcp_restart -i 3600 ckpt_*.dmtcp

# create output directory 

mkdir NLR_mapping_out_multihits

# create blast database 
makeblastdb -in EV_bedadeti.fna -dbtype nucl
makeblastdb -in EV_mazia.fna -dbtype nucl
makeblastdb -in Musa_acuminata.fna -dbtype nucl
makeblastdb -in Musa_balbisiana.fna -dbtype nucl

for i in $(cat genomes_db_list.txt)
 do
   for j in m$(cat NLRs_pred_list.txt)  # # NLRs_pred_list.txt: EV_mazia.NB-ARC.fasta, EV_beadeti.NB-ARC.fasta, musa_ac.NB-ARC.fasta, musa_ba.NB-ARC.fasta
     do
        # run tblastn 

        tblastn \
        -db ~/data/genomes/${i}.fna \
        -query ../pred_prot/PF00931_out/${j}.NB_ARC.fasta \ 
        -evalue 1e-5 \
        -outfmt 6  \
        -max_hsps 1 \
        -num_threads 48 \
        -out hmmNLR_mapping_out/${j}.${i}.e5.tblastn.out

       # extract hits

        # swap bed coordinates that start position > end
        awk '{if ($9 > $10) print $2"\t"$10"\t"$9"\t"$3; else print $2"\t"$9"\t"$10"\t"$3}' hmmNLR_mapping_out/${j}.${i}.e5.tblastn.out > hmmNLR_mapping_out/${j}.${i}.tblastn.0.bed

        # formatting bed file in case the new line character the end of each line
        cat hmmNLR_mapping_out/${j}.${i}.tblastn.0.bed  | awk '{print $1, $2, $3, $4}' OFS='\t' > hmmNLR_mapping_out/${j}.${i}.tblastn.1.bed 

        #  sort bed files
        bedtools sort -i hmmNLR_mapping_out/${j}.${i}.tblastn.1.bed > hmmNLR_mapping_out/${j}.${i}.tblastn.sorted.2.bed 

        # Merge coordinated
        bedtools merge -i hmmNLR_mapping_out/${j}.${i}.blastp.sorted.2.bed > hmmNLR_mapping_out/${j}.${i}.NLR_genes.tblasn.sorted.merged.bed
              
   done 
 done
 

# Extract non-overlapping coordinates of NLR loci between NLR encoding genes 
bedtools intersect -b NLR-annotator.out/mazia.NLR-annotator.out.bed -a hmmNLR_mapping_out/mazia.NLR_genes.tblasn.sorted.merged.bed -v > mazia_NLR_loci_non_overlapping.bed

# extract genomic sequence for blastx 
bedtools getfasta -fi EV_mazia.fna -bed mazia_NLR_loci_non_overlapping.bed > NhmmNLR_mapping_out/mazia.NLR_loci_non_overlapping.fasta

# EV_bdadeti 
bedtools intersect -a EV_bedadeti.fna -b bedadeti_NLR_loci_non_overlapping.bed > hmmNLR_mapping_out/bedadeti_non_overlapping.bed
bedtools getfasta -fi EV_bedadeti.fna -bed mazia_NLR_loci_non_overlapping.bed > hmmNLR_mapping_out/bedadeti.NLR_loci_non_overlapping.fasta

# Musa acuminata 
bedtools intersect -a Musa_balbisiana.fna -b musa_ac_NLR_loci_non_overlapping.bed > musa_ac_non_overlapping.bed

# Musa balbisiana 
bedtools intersect -a Musa_balbisiana.fna -b musa_ba_NLR_loci_non_overlapping.bed > musa_ba_non_overlapping.bed
bedtools getfasta -fi Musa_balbisiana.fna -bed musa_ba_NLR_loci_non_overlapping.bed > hmmNLR_mapping_out/musa_ba.NLR_loci_non_overlapping.fasta


# TBLASTX unmapped NLRs against NCBI non-redundant protein db

mkdir unmapped_NLR_blastx

for i in mazia bedadeti musa_ba
do
tblastx \
-db ~/NCBI_NT/BLASTDB/nt \
-query hmmNLR_mapping_out/musa_ba.NLR_loci_non_overlapping.fasta \
-evalue 1e-5 \
-outfmt 6  \
-max_hsps 1 \
-num_threads 16 \
-out unmapped_NLRs_blastx/${i}.unmapped_NLR.blastx.e5.out
done 

#!/bin/bash
#SBATCH --job-name=mac
#SBATCH --partition=cnode
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem-per-cpu=4571
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


# Extract Hidden Markov Models (HMM) of NLR from Pfam 

wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz 
tar xzfv Pfam-A.hmm.gz 
hmmpress Pfam-A.hmm

## Pfam IDs used for NLR: NB-ARC (PF00931), TIR (PF01582, PF13676), RPW8-type (PF05659), Late blight resistance protein R1 (PF12061), Rx-type (PF18052) and LRRs (PF00560, PF13516, PF18805, PF13855, PF07725, PF12799, PF01463, PF01462, PF18831, PF18837, PF08263, PF07723, and PF13306

hmmfetch -f Pfam_A/Pfam-A.hmm Pfam_NLR_ID.txt > Pfam_NLR.hmm 
hmmfetch -f Pfam_A/Pfam-A.hmm PF00931.txt > PF00931.hmm 

## NLR prediction  

#input_prots=$(cat prot_list.v1.txt )
input_prots=musa_ac
path_intpr=/home/lifesci/lfrwtp/apps/interproscan-5.52-86.0/

# create output directories 
mkdir interproscan_out AnnaRefPlant_NLR.out

# local diamond db of NLR from DOI: 10.1016/j.molp.2021.08.001 and https://doi.org/10.1371/journal.pbio.3001124

diamond makedb --in db_fasta/ANN_RePlant_NLR.fasta -d  db_fasta/ANN_RePlant_NLR.fasta

# predicted proteins of EV and Musa: EV_mazia.protein.fa, EV_bedadeti.protein.fa, musa_ac.protein.fa, musa_ba.protein.fa 
for i in EV_mazia EV_bedadeti musa_ac musa_ba
do
diamond blastp \
--query EV_Musa_proteins/${i}.protein.fa \
--db db_fasta/ANN_RePlant_NLR.fasta \
--ultra-sensitive \
--masking 0 \
--out AnnaRefPlant_NLR.out/${i}.blastp.out \
--outfmt 6  \
--compress 0  \
--evalue 1e-5 \
--threads 68

# extract hits

# swap bed coordinates that start position > end
awk '{if ($7 > $8) print $1"\t"$8"\t"$7"\t"$3; else print $1"\t"$7"\t"$8"\t"$3}' AnnaRefPlant_NLR.out/${i}.blastp.out > AnnaRefPlant_NLR.out/${i}.blastp.0.bed

# sort bed files
bedtools sort -i AnnaRefPlant_NLR.out/${i}.blastp.0.bed > AnnaRefPlant_NLR.out/${i}.blastp.1.sorted.bed

# formatting bed file in case the new line character the end of each line
cat AnnaRefPlant_NLR.out/${i}.blastp.1.sorted.bed | awk '{print $1, $2, $3, $4}' OFS='\t' > AnnaRefPlant_NLR.out/${i}.blastp.out.2.bed

#  sort bed files
bedtools sort -i AnnaRefPlant_NLR.out/${i}.blastp.out.2.bed > AnnaRefPlant_NLR.out/${i}.blastp.out.sorted.bed

# Merge coordinated
bedtools merge -i AnnaRefPlant_NLR.out/${i}.blastp.out.sorted.bed > AnnaRefPlant_NLR.out/${i}.blastp.sorted.merged.bed

#extract sequence
bedtools getfasta -fi EV_Musa_proteins/${i}.protein.faa -bed AnnaRefPlant_NLR.out/${i}.blastp.sorted.merged.bed > AnnaRefPlant_NLR.out/${i}.blastp.AnnaRefPlant_NLR.fasta


# generate alignment 
clustalo \
--in AnnaRefPlant_NLR.out/${i}.blastp.AnnaRefPlant_NLR.fasta \
--out AnnaRefPlant_NLR.out/${i}.blastp.AnnaRefPlant_NLR.aln.fasta \
--threads 28 \
--seqtype Protein


# hmm build
hmmbuild AnnaRefPlant_NLR.out/${i}.blastp.AnnaRefPlant_NLR.hmm AnnaRefPlant_NLR.out/${i}.blastp.AnnaRefPlant_NLR.aln.fasta

cat Pfam_NLR.hmm AnnaRefPlant_NLR.out/${i}.blastp.AnnaRefPlant_NLR.hmm > AnnaRefPlant_NLR.out/Pfam_RefPlant.${i}.hmm

# NLR prediction using combined HMM NLR domains 
## hmmsearch
hmmsearch --domtblout AnnaRefPlant_NLR.out/${i}.hmmbuild.hmmsearch_pfam_dbtblout_e03.out -E 0.0001 --domE 0.001 --cpu 28 -o AnnaRefPlant_NLR.out/${i}.hmmbuild.hmmsearch_pfam.e03.txt AnnaRefPlant_NLR.out/${i}.hmm EV_Musa_proteins/${i}.fasta

## extract envelope coordinate
grep -v '#' AnnaRefPlant_NLR.out/${i}.hmmbuild.hmmsearch_pfam_dbtblout_e03.out | awk '{print $1, $20, $21}' OFS='\t' > AnnaRefPlant_NLR.out/${i}.hmmbuild.hmmsearch_pfam_envelope_coordinates.bed

## sort bed files
bedtools sort -i AnnaRefPlant_NLR.out/${i}.hmmbuild.hmmsearch_pfam_envelope_coordinates.bed > AnnaRefPlant_NLR.out/${i}.hmmbuild.hmmsearch_pfam_envelope_coordinates.sorted.bed

## merge overlapping coordinates
bedtools merge -i AnnaRefPlant_NLR.out/${i}.hmmbuild.hmmsearch_pfam_envelope_coordinates.sorted.bed  > AnnaRefPlant_NLR.out/${i}.hmmbuild.hmmsearch_pfam_envelope_coordinates.sorted.merged.bed

## extract coordinates sequence
bedtools getfasta -fi EV_Musa_proteins/${i}.fasta -bed AnnaRefPlant_NLR.out/${i}.hmmbuild.hmmsearch_pfam_envelope_coordinates.sorted.merged.bed > AnnaRefPlant_NLR.out/${i}.hmmbuild.hmmsearch_pfam_envelope_coordinates.fasta


# Finetuning NLR prediction to those only containing NB-ARC domain 


## hmmsearch
hmmsearch --domtblout AnnaRefPlant_NLR.out/${i}.PF00931.hmmsearch_pfam_dbtblout_e03.out -E 0.0001 --domE 0.001 --cpu 28 -o AnnaRefPlant_NLR.out/${i}.PF00931.hmmsearch_pfam.e03.txt AnnaRefPlant_NLR.out/${i}.hmmbuild.hmmsearch_pfam_envelope_coordinates.fasta

## extract envelope coordinate
grep -v '#' AnnaRefPlant_NLR.out/${i}.PF00931.hmmsearch_pfam_dbtblout_e03.out | awk '{print $1, $20, $21}' OFS='\t' > AnnaRefPlant_NLR.out/${i}.PF00931.hmmsearch_pfam_envelope_coordinates.bed

## sort bed files
bedtools sort -i AnnaRefPlant_NLR.out/${i}.PF00931.hmmsearch_pfam_envelope_coordinates.bed > AnnaRefPlant_NLR.out/${i}.PF00931.hmmsearch_pfam_envelope_coordinates.sorted.bed

## merge overlapping coordinates
bedtools merge -i AnnaRefPlant_NLR.out/${i}.PF00931.hmmsearch_pfam_envelope_coordinates.sorted.bed  > AnnaRefPlant_NLR.out/${i}.PF00931.hmmsearch_pfam_envelope_coordinates.sorted.merged.bed

## extract coordinates sequence
bedtools getfasta -fi EV_Musa_proteins/${i}.fasta -bed AnnaRefPlant_NLR.out/${i}.hmmbuild.hmmsearch_pfam_envelope_coordinates.sorted.merged.bed > AnnaRefPlant_NLR.out/${i}.PF00931.hmmsearch_pfam_envelope_coordinates.fasta


## for predicted proteins

#mkdir interproscan_blastp_updated_hmmbuild
output_dir=interproscan_PF00931_hmmsearch

$path_intpr/interproscan.sh \
-i  AnnaRefPlant_NLR.out/${i}.PF00931.hmmsearch_pfam_envelope_coordinates.fasta \
-f gff3,TSV \
-d ${output_dir} \
-cpu  28

#rm -rf temp

# mast and fimo
fimo -o fimo.PF00931.${i}.out ~/apps/NLR-Annotator/bin/meme.xml AnnaRefPlant_NLR.out/${i}.PF00931.hmmsearch_pfam_envelope_coordinates.fasta
mast ~/apps/NLR-Annotator/bin/meme.xml AnnaRefPlant_NLR.out/${i}.PF00931.hmmsearch_pfam_envelope_coordinates.fasta -o mast.${i}.PF00931.out

# Extracting NB-ARC domaine using NCBI Conserved Domain Database search https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi


done



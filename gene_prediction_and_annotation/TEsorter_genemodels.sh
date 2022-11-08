#!/bin/bash
#SBATCH --job-name=tesort
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL,END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk


# TEsorter

# Predicted proteins of EV mazia and Bedadeti (AED<=0.25), MA and MB.

for i in EV_mazia EV_bedadeti Musa_acuminata Musa_balbisiana 
do 
TEsorter \
${i}.fasta \
--hmm rexdb-plant \
--seq-type  prot \
--processors 32 \
--pass2-rule 80-80-80 \
--prefix ${i}_TE.rexdb-plant
done





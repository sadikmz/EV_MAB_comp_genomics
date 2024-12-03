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
cpus=$cpus
BASENAME=$(echo $softmasked_genome | sed "s/.fna//g")

# Predicted proteins of EV landraces  [Mazia and Bedadeti (AED<=0.25)], Ensete glaucum (EG), Musa acuminata (MA) and Musa balbisiana (MB). Predicted proteins of Musa species and EG were obtained from https://banana-genome-hub.southgreen.fr/
# ls -lht ~/data | awk '{print $9}'
# EV_mazia.fna 
# EV_bedadeti.fna
# Musa_acuminata.fna 
# Musa_balbisiana.fna 

for i in ~/data/*.fna 
do 
	BASENAME=$(basename $i| sed 's/.fna//g')

	TEsorter \
	${i}.fna \
	--hmm rexdb-plant \
	--seq-type  prot \
	--processors 32 \
	--pass2-rule 80-80-80 \
	--prefix ${i}_TEsorter.rexdb-plant
done





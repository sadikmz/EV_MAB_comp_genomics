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





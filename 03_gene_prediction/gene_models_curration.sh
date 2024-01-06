#!/bin/bash
#!/bin/
#SBATCH --job-name=orthf
#SBATCH --partition=hmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL,END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk


set -eu

# OrthoFinder: Genomic sequence unmasked 

Rscript extract_longest_isoforms.R 

orthofinder -f EV_MAB -a 48

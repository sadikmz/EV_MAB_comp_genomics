#!/bin/bash
#SBATCH --job-name=phylo
#SBATCH --partition=hmem
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

## for ANN_RefPlant

clustalo \
--in EV_MAB.NB_ARC.RPW8_non_plantNBARC.fasta \
--out EV_MAB.NB_ARC.RPW8_non_plantNBARC.phy \
--outfmt=phylip \
--threads 48 \
--iterations 100 \
--seqtype=Protein


## Trim alignments 

~/apps/trimal-1.4.1/source/trimal -in EV_MAB.NB_ARC.RPW8_non_plantNBARC.phy -out EV_MAB.NB_ARC.RPW8_non_plantNBARC_gt_90.phy -phylip -gt 0.90 
#positions with gaps above 90% in aligned sequences were removed using 

# ~/apps/trimal-1.4.1/source/trimal -in EV_MAB.NB_ARC.RPW8_non_plantNBARC.phy -out EV_MAB.NB_ARC.RPW8_non_plantNBARC_gt_80.phy -phylip -gt 0.80 
# #positions with gaps above 20% in aligned sequences were removed using 

# ~/apps/trimal-1.4.1/source/trimal -in EV_MAB.NB_ARC.RPW8_non_plantNBARC.phy -out EV_MAB.NB_ARC.RPW8_non_plantNBARC_gt_50.phy -phylip -gt 0.50 
# #positions with gaps above 50% in aligned sequences were removed using 

# Run raxml 
raxmlHPC-PTHREADS-AVX2 -T 48 -s EV_MAB.NB_ARC.RPW8_non_plantNBARC_gt_90.phy -n EV_MAB.NB_ARC.RPW8_non_plantNBARC_90_RAxML -m PROTGAMMAAUTO -f a -# 1000 -x 123456 -p 246810 
#positions with gaps above 90% in aligned sequences were removed using 



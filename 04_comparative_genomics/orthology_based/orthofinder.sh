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

orthofinder -f EV_MAB -a 48 # EV_MAB contain longest predicted protein isoforms of EV, and Musa spp.

# Orthofinder was also tested in the following ways, but infered paralogue proteins appears to show high amino acid identity (>80%) with a considerable coverage (>50%).    

# # add predicted proteins monocot species downnloaded from uniprot 
# Ananas_comosus 
# Dioscorea_alata 
# Manihot_esculenta 
# Oryza_sativa 
# Sorghum_biocolor 
# Triticum_aestivum 
# Zea_mays


# Run orthofinder 
orthofinder -f EV_MAB_monocot -a 48 

# Extract inferted paralogous coding genes in EV and MAB, and combine EV paralogous into a single fasta file which total will be three fasta files: EV, MA and MB. 

# Re-run orthofinder with more stringent parameter options (including dimond blastp --ultra-sensitive)
orthofinder -d -t 48 -a 8 -M msa -A muscle -T raxml-ng -f orthofinder_filtered
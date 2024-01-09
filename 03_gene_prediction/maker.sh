#!/bin/bash
#SBATCH --job-name=maker
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL,END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk

module purge
module load GCC/11.3.0 OpenMPI/4.1.4
module load MAKER/3.01.04

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cpus=48
maxIntronLen=82825
MAKERDIR="mazia"
genotype=mazia

# round 1
maker -CTL 
# In the maker_opts.ctl file add absolute path to the following files and make the change indicated 
# genome=absolute_path_to_genome_file
# est=abslolute_path_to_trinity_pasa_RNA-seq_assembly 
# protein=absolute_path_to_protein_homology_fasta_file 
# cpu=1
# est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
# protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no

srun maker maker_exe.ctl maker_opts.ctl maker_bopts.ctl -qq -base ${MAKERDIR} -fix_nucleotid -cpus $cpus

gff3_merge round1_hints/${genotype}.est.*.gff -o round1_hints/genome.snap01.all.gff

# generate the necessary files for training SNAP:
cd snap.round1
maker2zff ../predicted.maker.gff3

# detect erroneous predictions
fathom genome.ann genome.dna -validate > snap_validate_output.txt

# Check errors
cat snap_validate_output.txt | grep "error" > snap.error.out
cat snap_validate_output.txt | grep 'errors' | awk '{print $1, $2}' >  snap.scaffolds_gene_moder_error.txt
cat snap_validate_output.txt | grep 'errors' | awk '{print $2}' > snap.gene_moder_error.txt

# remove
# remove the erroneous model from the file containing all gene models
#grep -vwE "MODEL24650" genome.ann > genome.ann01
grep -vwE -f snap.gene_moder_error.txt genome.ann > genome.ann1

# Validate if error exists
fathom genome.ann genome.dna -validate | grep "error" > snap.summary.txt
#cat round00.summary.txt #21005 genes, 21005 OK, 11017 warnings, 0 errors

#mkdir snap01
cd train_snap.round1

# create inputs for training SNAP
fathom ../genome.ann1 ../genome.dna -categorize 1000
fathom uni.ann uni.dna -export 1000 -plus
forge export.ann export.dna
hmm-assembler.pl ${genotype}.snap1 . > ${genotype}.snap1.hmm

### extract gff from mak.00.output
cat ../../predicted.maker.1.gff3 | grep -v '##FASTA' | grep '##gff-version 3\|^##\|^scf_6838' > ../../predicted.maker.1.noseq.gff3


# round2 

```bash
# In the maker_opts.ctl file add absolute path to the following files and make the change indicated 
# genome=absolute_path_to_genome_file
# est=abslolute_path_to_trinity_pasa_RNA-seq_assembly 
# protein=absolute_path_to_protein_homology_fasta_file 
# cpu=1
# est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
# protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
# maker_gff=absolute_path_to_predicted.maker.1.noseq.gff3 #MAKER derived GFF3 file
# snaphmm=${genotype}.snap1.hmm #SNAP HMM file

srun maker maker_exe.ctl maker_opts.ctl maker_bopts.ctl -qq -base ${MAKERDIR} -fix_nucleotid -cpus $cpus

gff3_merge round1_hints/${genotype}.est.*.gff -o round1_hints/genome.snap01.all.gff

# train SNAP as described above and this will generate: {genotype}.snap2.hmm 

# round3: Repeat round 2 with new {genotype}.snap2.hmm. This was run for one more rounds so that snap trained 3x.
```

# round5
```
# In the maker_opts.ctl file add absolute path to the following files and make the change indicated 
# genome=absolute_path_to_genome_file
# est=abslolute_path_to_trinity_pasa_RNA-seq_assembly 
# protein=absolute_path_to_protein_homology_fasta_file 
# cpu=1
# est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
# protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
# maker_gff=absolute_path_to_predicted.maker.1.noseq.gff3 #MAKER derived GFF3 file
# snaphmm=${genotype}.snap3.hmm #SNAP HMM file
# gmhmm=path_to_braker_gmhmm.mod #GeneMark HMM file
# augustus_species=path_to_ensete_ventricosum_augustus #Augustus gene prediction species model

srun maker maker_exe.ctl maker_opts.ctl maker_bopts.ctl -qq -base ${MAKERDIR} -fix_nucleotid -cpus $cpus
```

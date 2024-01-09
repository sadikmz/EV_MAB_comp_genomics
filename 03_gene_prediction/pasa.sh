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

module load intel impi imkl
module purge

genotype=mazia 
trinity_rnaseq=data/Trinity.fasta
trinity_cleaned=data/Trinity.cleaned.fasta 
trans_gtf=${genotype}_rnaseq.gff
genome=${genotype}.fna 
cpus=48


~/apps/PASApipeline-v2.5.1/bin/seqclean /data/Trinity.fasta -c 48 -o $trinity_rnaseq -n 10000

# landrace/genotype mazia and bedadeti
~/apps/PASApipeline-v2.5.1/Launch_PASA_pipeline.pl \
-c ~/apps/PASApipeline-v2.5.1//pasa_conf/pasa.config \
-C -R \
-g ${genotype}.fna \
-t $trinity_cleaned \
--trans_gtf $trans_gtf \
--TRANSDECODER \
--ALIGNERS blat,gmap \
--CPU $cpus 


# pasa config detail 

## templated variables to be replaced exist as <__var_name__>

# database settings
#DATABASE=/tmp/hsperfdata_lfrwtp
DATABASE=/tmp/pasa_genotype
#DATABASE=/scratch/lifesci/lfrwtp/pasa_tmp
#######################################################
# Parameters to specify to specific scripts in pipeline
# create a key = "script_name" + ":" + "parameter" 
# assign a value as done above.

#spliced multiple mappings allowed
run_spliced_aligners.pl:-N=5

#script validate_alignments_in_db.dbi
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=80
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=80
validate_alignments_in_db.dbi:--MAX_INTRON_LENGTH=20000

#script subcluster_builder.dbi
subcluster_builder.dbi:-m=50

#annotation comparison
cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP=70
cDNA_annotation_comparer.dbi:--MIN_PERCENT_PROT_CODING=40
cDNA_annotation_comparer.dbi:--MIN_PERID_PROT_COMPARE=70
cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_FL_COMPARE=50
cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_NONFL_COMPARE=50
cDNA_annotation_comparer.dbi:--MIN_PERCENT_ALIGN_LENGTH=50
cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP_GENE_REPLACE=80
cDNA_annotation_comparer.dbi:--MAX_UTR_EXON


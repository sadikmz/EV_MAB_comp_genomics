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

module load GCC/8.3.0 OpenMPI/3.1.4

#~/apps/PASApipeline-v2.5.1/bin/seqclean /gpfs/home/lifesci/lfrwtp/analysis/annotation/mazia/Trinity.fasta -c 48 -o Trinity.cleaned.fasta -n 10000

#/home/lifesci/lfrwtp/miniconda3/envs/pasa/opt/pasa-2.3.3/scripts/create_sqlite_cdnaassembly_db.dbi \
#-c /home/lifesci/lfrwtp/miniconda3/envs/pasa/opt/pasa-2.3.3/pasa_conf/pasa.mazia.config \
#-S '/home/lifesci/lfrwtp/miniconda3/envs/pasa/opt/pasa-2.3.3/schema/cdna_alignment_sqliteschema
#
#/home/lifesci/lfrwtp/miniconda3/envs/pasa/opt/pasa-2.3.3/scripts/upload_transcript_data.dbi \
#-M '/home/lifesci/lfrwtp/analysis/annotation/mazia/pasa/tmp/maz.hints_mydb_pasa.sqlite' \
#-t data/Trinity.cleaned.fasta \
#-f NULL

#/home/lifesci/lfrwtp/apps/PASApipeline-v2.5.1/Launch_PASA_pipeline.pl \
/home/lifesci/lfrwtp/miniconda3/envs/funannotate/opt/pasa-2.4.1/Launch_PASA_pipeline.pl \
-c /home/lifesci/lfrwtp/miniconda3/envs/pasa/opt/pasa-2.3.3/pasa_conf/pasa.mazia.config \
-C -R \
-g ${mazia}.fna \
-t data/Trinity.cleaned.fasta \
--trans_gtf mazia_rnaseq.1.gff  \
--TRANSDECODER \
--ALIGNERS blat,gmap \
--CPU 48 
#--INVALIDATE_SINGLE_EXON_ESTS 


MAKERDIR="mazia"

~/apps/maker/bin/maker -cpus 48 -qq -base ${MAKERDIR} -fix_nucleotides


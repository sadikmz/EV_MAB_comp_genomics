#!/bin/bash

# when error with PERL5LIB user PerlLib from PASA
# PERL5LIB="$/home/lifesci/lfrwtp/apps/PASApipeline.v2.4.1/PerlLib"
# PERL5LIB="$/home/u1866313/apps/PASApipeline.v2.4.1/PerlLib"

## PASApipeline
/home/u1866313/apps/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl \
--ALIGNERS gmap,blat \
--MAX_INTRON_LENGTH 100000 \
--IMPORT_CUSTOM_ALIGNMENTS_GFF3 \
--trans_gtf cufflinks or stringtie--generated transcripts \
--run \
--annot_compare \
--ALT_SPLICE \
--genome \
--transcripts \
--TDN \
-T \
--TRANSDECODER \
--CPU \
--PASACONF $PASAHOME/pasa_conf/conf.txt \


/home/u1866313/apps/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl \
-c /home/u1866313/apps/PASApipeline.v2.4.1/pasa_conf/pasa.alignAssembly.conf \
-C -R \
-g $flye_v1_1 \
-t Trinity.fasta \
--ALIGNERS blat,gmap \
--CPU 32
perl Makefile.PL --mysql_config=/usr/local/bin/mysql_config

/Users/u1866313/Apps/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl \
-c alignAssembly.config \
-C -R \
-g genome_sample.fasta \
-t Trinity.fasta \
--ALIGNERS blat,gmap

## configure zlib

export PATH=$PATH:$home/u1866313/miniconda3/envs/zlib/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$home/u1866313/miniconda3/envs/zlib/lib/
export LIBRARY_PATH=$LIBRARY_PATH:$HOME/Softwares/library/Zlib/zlib-1.2.11/lib/
export C_INCLUDE_PATH=$home/u1866313/miniconda3/envs/zlib/include/
export CPLUS_INCLUDE_PATH=$home/u1866313/miniconda3/envs/zlib/include/
export PKG_CONFIG_PATH=$home/u1866313/miniconda3/envs/zlib/lib/pkgconfig/


./__run_sample_pipeline.pl \
--align_assembly_config sqlite.confs/alignAssembly.config \
--annot_compare_config sqlite.confs/annotCompare.config  \
--TRANSDECODER

(pasa0) [lfrwtp@orac:login1 pasa_conf]$ cat ../sample_data0000/sqlite.confs/alignAssembly.config

## templated variables to be replaced exist as <__var_name__>

# Pathname of an SQLite database
# If the environment variable DSN_DRIVER=mysql then it is the name of a MySQL database
DATABASE=/tmp/sample_mydb_pasa.sqlite


#######################################################
# Parameters to specify to specific scripts in pipeline
# create a key = "script_name" + ":" + "parameter"
# assign a value as done above.

#script validate_alignments_in_db.dbi
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=75
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=95
validate_alignments_in_db.dbi:--NUM_BP_PERFECT_SPLICE_BOUNDARY=0

#script subcluster_builder.dbi
subcluster_builder.dbi:-m=50

(pasa0) [lfrwtp@orac:login1 pasa_conf]$ cat ../sample_data0000/sqlite.confs/annotCompare.config

## templated variables to be replaced exist as <__var_name__>

# Pathname of an SQLite database
# If the environment variable DSN_DRIVER=mysql then it is the name of a MySQL database
DATABASE=/tmp/sample_mydb_pasa.sqlite
DATABASE=/tmp/sample_mydb_pasa.sqlite


#######################################################
# Parameters to specify to specific scripts in pipeline
# create a key = "script_name" + ":" + "parameter"
# assign a value as done above.


#script cDNA_annotation_comparer.dbi
cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP=<__MIN_PERCENT_OVERLAP__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_PROT_CODING=<__MIN_PERCENT_PROT_CODING__>
cDNA_annotation_comparer.dbi:--MIN_PERID_PROT_COMPARE=<__MIN_PERID_PROT_COMPARE__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_FL_COMPARE=<__MIN_PERCENT_LENGTH_FL_COMPARE__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_NONFL_COMPARE=<__MIN_PERCENT_LENGTH_NONFL_COMPARE__>
cDNA_annotation_comparer.dbi:--MIN_FL_ORF_SIZE=<__MIN_FL_ORF_SIZE__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_ALIGN_LENGTH=<__MIN_PERCENT_ALIGN_LENGTH__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP_GENE_REPLACE=<__MIN_PERCENT_OVERLAP_GENE_REPLACE__>
cDNA_annotation_comparer.dbi:--STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE=<__STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE__>
cDNA_annotation_comparer.dbi:--TRUST_FL_STATUS=<__TRUST_FL_STATUS__>
cDNA_annotation_comparer.dbi:--MAX_UTR_EXONS=<__MAX_UTR_EXONS__>
cDNA_annotation_comparer.dbi:--GENETIC_CODE=<__GENETIC_CODE__>

#/home/lifesci/lfrwtp/apps/PASApipeline.v2.4.1/pasa_conf/annotCompare.config
/home/lifesci/lfrwtp/apps/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl \
-c alignAssembly.config \
-C -R \
-g /home/lifesci/lfrwtp/data/annotation/rnaseq/epo.flye.v1.scaffolds.softmaked.fna \
-t /home/lifesci/lfrwtp/data/annotation/rnaseq/Trinity.fasta \
--ALIGNERS blat,gmap \
--CPU 68 \

cufflinks --total-hits-norm  --num-threads 16 --output-dir epo_mazia_cufflinks




## De novo assembly transcriptom using Trinity 

# RNA-seq assembly using Trinity 
```bash
#!/bin/bash

# get Trinity assembly from https://github.com/trinityrnaseq/trinityrnaseq/releases

export TRINITY_HOME=/home/u1866313/apps/trinityrnaseq-v2.15.1

RNAseq_DIR=path_to_rnaseq_reads_dir

# get trim_galore from https://github.com/FelixKrueger/TrimGalore

trim_galore -cores 48 --quality 30 --output_dir trim_out --paired $RNAseq_DIR/DP8400009186BL_L01_573_?.fq.gz
trim_galore -cores 48 --quality 30 --output_dir trim_out --paired $RNAseq_DIR/DP8400009186BL_L01_574_?.fq.gz
trim_galore -cores 48 --quality 30 --output_dir trim_out --paired $RNAseq_DIR/DP8400009186BL_L01_575_?.fq.gz


export PERL5LIB=""

## Run Trinity 
$TRINITY_HOME/Trinity \
--seqType fq \
--left ${RNAseq_DIR}/DP8400009186BL_L01_573_1_val_1.fq.gz,${RNAseq_DIR1}/DP8400009186BL_L01_574_1_val_1.fq.gz,${RNAseq_DIR1}/DP8400009186BL_L01_575_1_val_1.fq.gz \
--right ${RNAseq_DIR}/DP8400009186BL_L01_573_2_val_2.fq.gz,${RNAseq_DIR1}/DP8400009186BL_L01_574_2_val_2.fq.gz,${RNAseq_DIR1}/DP8400009186BL_L01_575_2_val_2.fq.gz \
--SS_lib_type RF  \
--min_contig_length 100 \
--max_memory 300G \
--CPU 36 \
--full_cleanup \
--output trinity_out 
```
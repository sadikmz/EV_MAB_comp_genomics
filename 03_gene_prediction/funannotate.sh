funannotate predict \
--input data/epo.flye.v1.scaffolds.softmaked.fna \
--out funannotate.out \
--species "Ensete ventricosum" -a Ensete_ventricosum.parameters.json \
--transcript_evidence data/Trinity.fasta \
--rna_bam data/flye.v1.1.rnaseq.sorted.bam \
--pasa_gff data/mazia_trinity_epo.hints_mydb_pasa.sqlite.pasa_assemblies.gff3:10 \
--other_gff data/final_annotation.gff:7 \
--organism other \
--busco_db embryophyta \
--augustus_gff  data/brakerv2.16.augustus_addUTR.hints_utr.gtf:7 \
--genemark_gtf data/genemark.d.gtf:5 \
--cpus 28 \
--max_intronlen 19000 \
--weights snap:5


# train

funannotate train \
--input data/epo.flye.v1.scaffolds.softmaked.fna \
--out funannot.out \
--left data/rnaseq/DP8400009186BL_L01_573_1.fq.gz data/rnaseq/DP8400009186BL_L01_574_1.fq.gz data/rnaseq/DP8400009186BL_L01_575_1.fq.gz \
--right data/rnaseq/DP8400009186BL_L01_573_2.fq.gz data/rnaseq/DP8400009186BL_L01_574_2.fq.gz data/rnaseq/DP8400009186BL_L01_575_2.fq.gz \
--no_trimmomatic \
--trinity data/Trinity.fasta \
--species "Ensete ventricosum" \
--cpus 28


funannotate predict \
--input data/epo.flye.v1.scaffolds.softmaked.fna \
--out funannotate.out \
--cpus 28 \
--augustus_species ensete_ventricosum \
--transcript_evidence data/Trinity.fasta \
--rna_bam data/flye.v1.1.rnaseq.sorted.bam \
--pasa_gff pasa/Ensete_ventricosum_pasa.assemblies.fasta.transdecoder.genome.gff3:10 \
--other_gff funannotate/gemoma.gff3:7 \
--organism other \
--ploidy 2 \
--repeats2evm \
--busco_db embryophyta \
--repeat_filter  overlap blast \
--genemark_gtf data/genemark.d.gtf --max_intronlen 19000 \
--weights snap:1 genemark:2 augustus:7 glimmerhmm:1

funannotate predict \
--input data/epo.flye.v1.scaffolds.softmaked.fna \
--out funannotate.out \
--cpus 28 \
--augustus_species ensete_ventricosum \
--protein_evidence data/enseteg.musa.ensemble_mono.prot.fa \
--transcript_evidence data/Trinity.fasta \
--rna_bam data/flye.v1.1.rnaseq.sorted.bam \
--pasa_gff funannotate_train.pasa.gff3:10 \
--other_gff funannotate/gemoma.gff3:7 \
--organism other \
--ploidy 2 \
--repeats2evm \
--busco_db embryophyta \
--repeat_filter  overlap blast \
--genemark_gtf data/genemark.d.gtf --max_intronlen 19000 \
--weights snap:1 genemark:4 augustus:7 glimmerhmm:1

#!/bin/bash

# Gemoma with protein coding genes 

## Installation 

conda create -n gemoma -c bioconda gemoma

~/miniconda3/envs/gemoma/share/gemoma-1.7.1-0/GeMoMa-1.7.1.jar

## Download assembly and gff files of species closely related to the your genome from Phytozome v13 in directory named gemoma_list
## Things to modity in the gemoma script 
### "i" is just abbrivation of the species name and I took the first and letter of species and genus name for simplicity. Modify your own accordingly  
### path to input genome in "genome" and "g", and "a" for path to annotation gff file


## path to repeatmasked genome to be annotated 

genome=../repeatmasked_genome.fna 

# create output directory 

dir=GeMoMa_annotation.out

# activate Gemoma's installation conda environment 

conda activate gemoma 

# Run GeMoMA

java -jar -Xms10G -Xmx100G ~/miniconda3/envs/gemoma/share/gemoma-1.7.1-0/GeMoMa-1.7.1.jar CLI GeMoMaPipeline \
t=$genome \
d=DENOISE AnnotationFinalizer.r=NO o=true GeMoMa.s=false GeMoMa.t=3600 sc=true p=true pc=true tblastn=true \
s=own i=MA g=../gemoma_list/Musa_acuminata_pahang_v4.fna a=../gemoma_list/Musa_acuminata_pahang_v4.gff3 \
s=own i=MB g=../gemoma_list//Musba_assembly.fna a=../gemoma_list/Musba_assembly.gff3 \
s=own i=MC g=../gemoma_list/banksii_final_assembly_v2_filtered.fasta a=../gemoma_list/Maban_annoted.gff3 \
s=own i=EG g=../gemoma_list/ensete_glaucum.assembly.fna a=../gemoma_list/ensete_glaucum.assembly.fna.gff3 \
s=own i=AT g=../gemoma_list/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa a=../gemoma_list/Arabidopsis_thaliana.TAIR10.51.gff3 \
s=own i=AC g=../gemoma_list/Ananas_comosus.F153.dna.toplevel.fa a=../gemoma_list/Ananas_comosus.F153.51.gff3 \
s=own i=BD g=../gemoma_list/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.dna_sm.toplevel.fa a=../gemoma_list/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.51.gff3 \
s=own i=DR g=../gemoma_list/Dioscorea_rotundata.TDr96_F1_v2_PseudoChromosome.dna_sm.toplevel.fa a=../gemoma_list/Dioscorea_rotundata.TDr96_F1_v2_PseudoChromosome.51.gff3 \
s=own i=ME g=../gemoma_list/Manihot_esculenta.Manihot_esculenta_v6.dna_sm.toplevel.fa a=../gemoma_list/Manihot_esculenta.Manihot_esculenta_v6.51.gff3 \
s=own i=OI g=../gemoma_list/Oryza_indica.ASM465v1.dna_sm.toplevel.fa a=../gemoma_list/Oryza_indica.ASM465v1.51.sortd.gff3 \
s=own i=OS g=../gemoma_list/Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa a=../gemoma_list/Oryza_sativa.IRGSP-1.0.51.gff3 \
s=own i=SB g=../gemoma_list/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna_sm.toplevel.fa a=../gemoma_list/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.51.gff3 \
s=own i=ZM g=../gemoma_list/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna_sm.toplevel.fa a=../gemoma_list/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.51.gff3 \
# restart=true \
GeMoMa.prefix=gemoma \
threads=28 \
outdir=${dir} 

# filter prediction 

cat $dir/unfiltered_predictions_from_species_?.gff >  unfiltered_predictions_from_species_all.gff

java -jar -Xms10G -Xmx250G ../GeMoMa/GeMoMa-1.7.1.jar CLI GAF \
p=pred \
g=unfiltered_predictions_from_all_species.gff \
outdir=GeMoMa_filtered_annotation.out 


# Gemoma with RNA-seq and protein coding genes
GeMoMa GeMoMaPipeline \
t=target_genome \
q=query proteins \
r=mapped \
ERE.m=RNASeq.sorted.bam \
ERE.c=true \
tblastn \
d=DENOISE \ # remove questionable interons
GeMoMa.s \ #
AnnotationFinalizer.u \
AnnotationFinalizer.r \
GAF.c \
GAF.f \
GAF.a \
threads 64 \
GeMoMa.t 7200 \
s=own i=AT g=gemoma_list/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa a=gemoma_list/Arabidopsis_thaliana.TAIR10.51.gff3 \
s=own i=AC g=gemoma_list/Ananas_comosus.F153.dna.toplevel.fa a=gemoma_list/Ananas_comosus.F153.51.gff3 \
s=own i=BD g=gemoma_list/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.dna_sm.toplevel.fa a=gemoma_list/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.51.gff3
s=own i=DR g=Dioscorea_rotundata.TDr96_F1_v2_PseudoChromosome.dna_sm.toplevel.fa a=Dioscorea_rotundata.TDr96_F1_v2_PseudoChromosome.51.gff3
s=own i=ME g=Manihot_esculenta.Manihot_esculenta_v6.dna_sm.toplevel.fa a=Manihot_esculenta.Manihot_esculenta_v6.51.gff3
s=own i=OI g=Oryza_indica.ASM465v1.dna_sm.toplevel.fa a=Oryza_indica.ASM465v1.51.sortd.gff3
s=own i=OS g=Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa a=Oryza_sativa.IRGSP-1.0.51.gff3
s=own i=SB g=Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna_sm.toplevel.fa a=Sorghum_bicolor.Sorghum_bicolor_NCBIv3.51.gff3
s=own i=ZM g=Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna_sm.toplevel.fa a=Zea_mays.Zm-B73-REFERENCE-NAM-5.0.51.gff3 \
outdir GeMoMA.out


java -jar -Xms10G -Xmx250G ../GeMoMa/GeMoMa-1.7.1.jar CLI GAF \
p=pred g=unfiltered_predictions_from_species_0.gff,unfiltered_predictions_from_species_1.gff,unfiltered_predictions_from_species_2.gff,unfiltered_predictions_from_species_3.gff,unfiltered_predictions_from_species_4.gff,unfiltered_predictions_from_species_5.gff,unfiltered_predictions_from_species_6.gff,unfiltered_predictions_from_species_7.gff,unfiltered_predictions_from_species_8.gff \
outdir=GeMoMA_fileter.out \

# UTR prediction
java -jar -Xms10G -Xmx250G ../GeMoMa/GeMoMa-1.7.1.jar CLI AnnotationFinalizer \
g=$genome \
a=annotation.gff \
u=YES \
c=UNSTRANDED \
coverage_unstranded=../GeMoMa_temp/GeMoMaPipeline-18391413622443481850/coverage_2.bedgraph \
i=/home/u1866313/Novogene/annotation/flye.v1/homology_hints/GeMoMa_temp/GeMoMaPipeline-18391413622443481850/introns.gff \
outdir=GeMoMA_fileter_addUTR.out
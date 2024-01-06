load("bedadeti_EG_MABS_alignment_coverage_dedup.RData")

# bedadti_musa_ac_align_cov <- mazia_musa_ac_align_cov 
# bedadeti_musa_ba_align_cov <- mazia_musa_ba_align_cov
bedadeti_musa_sc_align_cov <- mazia_musa_sc_align_cov
# bedadti_musa_ac_align_cov <- mazia_ensete_gl_align_cov

rm (mazia_musa_ac_align_cov,mazia_musa_sc_align_cov,mazia_ensete_gl_align_cov, mazia_musa_ba_align_cov)

bedadeti_musa_sc_align_cov %>%
  head()

## format input data and generate grange 

# bedadeti 

bedadti_musa_sc_align_cov_grange <-
  bedadeti_musa_sc_align_cov %>% 
  # dplyr::rename(y = coverage) %>% 
  filter(str_detect(chr,"^chr")) %>%
  filter(coverage > 0) %>%
  mutate(
    # chr = str_replace(chr, "Bchr", "chr"),
    start = as.numeric(start),
    coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
                         TRUE ~ as.numeric(round(log(coverage),2)))) %>%
  bed_to_granges ()

bedadti_musa_sc_align_cov_grange %>%
  head()

rm(bedadeti_musa_sc_align_cov)

# Mazia

load("mazia_EG_MABS_alignment_coverage_dedup.RData")

rm(mazia_ensete_gl_align_cov,mazia_musa_ac_align_cov,mazia_musa_ba_align_cov)

# Generate grange 

mazia_musa_sc_align_cov_filtered_grange <-
  mazia_musa_sc_align_cov %>% 
  # dplyr::rename(y = coverage) %>%
  filter(str_detect(chr,"^chr")) %>%
  filter(coverage > 0) %>%
  mutate(
    # chr = str_replace(chr, "Bchr", "chr"),
    start = as.numeric(start),
    coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
                         TRUE ~ as.numeric(round(log(coverage),2)))) %>%
  bed_to_granges ()



mazia_musa_sc_align_cov_filtered_grange %>%
  tail

rm (mazia_musa_sc_align_cov)

# Generate karyotype plot

library(karyoploteR)
library(GenomicRanges)


MS_genome_bed <- 
  read.delim("Musa_schizocarpa.genome.bed") %>%
  filter(str_detect(seq.name,"chr")) %>%
  dplyr::rename(chr = seq.name)

centro_density <- 
  
  read.delim("ALL_centro.density.100kb", header= F, sep = ' ') %>%
  dplyr::rename(chr = V1, 
                start = V2, 
                end = V3, 
                density = V4) %>%  
  dplyr::filter(str_detect(chr,"^S")) %>%
  mutate(chr = str_replace(chr, "S","chr")) %>% 
  bed_to_granges () 


kp <- plotKaryotype(genome=MS_genome_bed, ideogram.plotter = NULL, plot.type=6 )
# kpAddCytobandsAsLine(kp)
kpBars(kp, mazia_musa_sc_align_cov_filtered_grange,
       y1=mazia_musa_sc_align_cov_filtered_grange$id,
       ymax=max(mazia_musa_sc_align_cov_filtered_grange$id)/2,
       col="#006400", r0=0.55, r1=1, border=NA)

kpBars(kp, centro_density,
       y1=centro_density$id,
       ymax=max(centro_density$id)/0.2,
       col="red", r0=0.42, r1=1,border=NA)

kpBars(kp, bedadti_musa_sc_align_cov_grange,
       y1=bedadti_musa_sc_align_cov_grange$id,
       ymax=max(bedadti_musa_sc_align_cov_grange$id)/3.2,
       col="#0000FF", r0=0.35, r1=0,border=NA)


library(karyoploteR)

gff.file <- "http://plasmodb.org/common/downloads/release-37/PvivaxP01/gff/data/PlasmoDB-37_PvivaxP01.gff"

header.lines <- readLines(gff.file, n = 30)

#The lines with the standard chromosomes start with "##sequence-region PvP01".
#Select them.
ll <- header.lines[grepl(header.lines, pattern = "##sequence-region PvP01")]

#split them by space, and create a data.frame
gg <- data.frame(do.call(rbind, strsplit(ll, split = " ")))
gg[,3] <- as.numeric(as.character(gg[,3]))
gg[,4] <- as.numeric(as.character(gg[,4]))

#and create a GRanges with the information
PvP01.genome <- toGRanges(gg[,c(2,3,4)])

PvP01.genome


kp <- plotKaryotype(genome=PvP01.genome)

kp <- plotKaryotype(genome=PvP01.genome, ideogram.plotter = NULL)
kpAddCytobandsAsLine(kp)



library(rtracklayer)
features <- import(gff.file)
table(features$type)

genes <- features[features$type=="gene"]


kp <- plotKaryotype(genome=PvP01.genome, ideogram.plotter = NULL)
kpAddCytobandsAsLine(kp)
kpPlotRegions(kp, data=genes)


kp <- plotKaryotype(genome=PvP01.genome, ideogram.plotter = NULL, plot.type=2)
kpAddCytobandsAsLine(kp)
kpPlotRegions(kp, data=genes[strand(genes)=="+"], avoid.overlapping = FALSE)
kpPlotRegions(kp, data=genes[strand(genes)=="-"], avoid.overlapping = FALSE, data.panel=2)



pp <- getDefaultPlotParams(plot.type=2)
pp$data1outmargin <- 100
pp$data2outmargin <- 100
pp$topmargin <- 450
kp <- plotKaryotype(genome=PvP01.genome, ideogram.plotter = NULL, plot.type=2, plot.params = pp)
kpAddCytobandsAsLine(kp)
kpAddMainTitle(kp, "Plasmodium Vivax - PvP01 with genes", cex=2)
kpPlotRegions(kp, data=genes[strand(genes)=="+"], avoid.overlapping = FALSE, col="deepskyblue")
kpPlotRegions(kp, data=genes[strand(genes)=="-"], avoid.overlapping = FALSE, col="gold", data.panel=2)
kpAddLabels(kp, "strand +", cex=0.8, col="#888888")
kpAddLabels(kp, "strand -", data.panel=2, cex=0.8, col="#888888")



# MA tranalated genome NLR

## Genome bed coordinate 

MA_genome_bed  <-
  read.delim("~/rstudio/bgimazia_musa/NLR_translated_genome/karyotype/Musa_acuminata.genome.bed", header = T) %>%
  dplyr::rename(chr = seq.name) %>%
  filter(!str_detect(chr, "putative_mitochondrion")) %>%
  mutate(end = str_replace(end, "41765374",  "48736620"))

## NLR  coordinate 


MA_predictedProteins_NLR_bed <- read.delim("~/rstudio/bgimazia_musa/NLR_predicted_proteins/karyotype/Musa_acuminata.tblastn.sorted.merged.bed", header = F)  %>%
  dplyr::rename(chr = V1,
                start  = V2,  
                end = V3)  %>%
  mutate(chr = str_remove(chr, "B"))


MA_predictedProteins_NLR_bed %>% head()
MA_translated_NLR_bed <- read.delim("~/rstudio/bgimazia_musa/NLR_translated_genome/karyotype/Musa_acuminata.tblastn.sorted.merged.bed", header = F)  %>%
  dplyr::rename(chr = V1,
                start  = V2,  
                end = V3)
MA_translated_NLR_bed %>% head()


#Open the device
tiff("result_plots/MA_predicted_vs_translated_NLR_karyoplot.tiff", width = 800, height = 400)

# Generate karyoplot 
pp <- getDefaultPlotParams(plot.type=2)
pp$data1outmargin <- 10
pp$data2outmargin <- 10
pp$topmargin <- 400
kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=2, plot.params = pp, cex= 1.5)
kpAddCytobandsAsLine(kp)
kpAddMainTitle(kp, "M. acuminata", cex=2)
kpPlotRegions(kp, data=MA_predictedProteins_NLR_bed, avoid.overlapping = FALSE, col="deepskyblue")
kpPlotRegions(kp, data=MA_translated_NLR_bed, avoid.overlapping = FALSE, col="red", data.panel=2)
# kpAddLabels(kp, "MA", cex=0.8, col="#888888")
# kpAddLabels(kp, "MB", data.panel=2, cex=0.8, col="#888888")


#Close the device
dev.off()

# MA annotated proteins NLR

## Genomic coordinate
MB_genome_bed  <-
  read.delim("~/rstudio/bgimazia_musa/NLR_translated_genome/karyotype/Musa_balbisiana.genome.bed", header = T) %>%
  dplyr::rename(chr = seq.name) %>%
  filter(!str_detect(chr, "Bscaffold")) %>%
  mutate(chr = str_remove(chr, "B"))

# NLRs

MB_translated_NLR_bed <- read.delim("~/rstudio/bgimazia_musa/NLR_translated_genome/karyotype/Musa_balbisiana.tblastn.sorted.merged.bed", header = F)  %>%
  dplyr::rename(chr = V1,
                start  = V2,  
                end = V3)  %>%
  mutate(chr = str_remove(chr, "B"))

MB_predictedProteins_NLR_bed <- read.delim("~/rstudio/bgimazia_musa/NLR_predicted_proteins/karyotype/Musa_balbisiana.tblastn.sorted.merged.bed", header = F)  %>%
  dplyr::rename(chr = V1,
                start  = V2,  
                end = V3)  %>%
  mutate(chr = str_remove(chr, "B"))


#Open the device
tiff("result_plots/MB_predicted_vs_translated_NLR_karyoplot.tiff", width = 800, height = 400)


# Generate karyoplot 
pp <- getDefaultPlotParams(plot.type=2)
pp$data1outmargin <- 100
pp$data2outmargin <- 100
pp$topmargin <- 450
kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=2, plot.params = pp, cex=1.5)
kpAddCytobandsAsLine(kp)
kpAddMainTitle(kp, "M. balbisiana", cex=2)
kpPlotRegions(kp, data=MB_predictedProteins_NLR_bed, avoid.overlapping = FALSE, col="deepskyblue")
kpPlotRegions(kp, data=MB_translated_NLR_bed, avoid.overlapping = FALSE, col="red", data.panel=2)
# kpAddLabels(kp, "MA", cex=0.8, col="#888888")
# kpAddLabels(kp, "MB", data.panel=2, cex=0.8, col="#888888")


#Close the device
dev.off()


# MS
# Genome bed coordinate 

MS_genome_bed  <-
  read.delim("~/rstudio/bgimazia_musa/NLR_translated_genome/karyotype/Musa_schizocarpa.genome.bed", header = T) %>%
  dplyr::rename(chr = seq.name) %>%
  filter(!str_detect(chr, "Mschizocarpa_scaffold"),
         !str_detect(chr,"chloro"),
         !str_detect(chr,"mito")) 
# NLRs
MS_translated_NLR_bed <- read.delim("~/rstudio/bgimazia_musa/NLR_translated_genome/karyotype/Musa_shizocarpa.tblastn.sorted.merged.bed", header = F)  %>%
  dplyr::rename(chr = V1,
                start  = V2,  
                end = V3)

MS_predictedProteins_NLR_bed <- read.delim("~/rstudio/bgimazia_musa/NLR_predicted_proteins/karyotype/Musa_shizocarpa.tblastn.sorted.merged.bed", header = F)  %>%
  dplyr::rename(chr = V1,
                start  = V2,  
                end = V3)

## 
#Open the device
tiff("result_plots/MS_predicted_vs_translated_NLR_karyoplot.tiff", width = 800, height = 400)

pp <- getDefaultPlotParams(plot.type=2)
pp$data1outmargin <- 100
pp$data2outmargin <- 100
pp$topmargin <- 450
kp <- plotKaryotype(genome=MS_genome_bed, ideogram.plotter = NULL, plot.type=2, plot.params = pp, cex = 1.5)
kpAddCytobandsAsLine(kp)
kpAddMainTitle(kp, "M. schizocarpa", cex=2)
kpPlotRegions(kp, data=MS_predictedProteins_NLR_bed, avoid.overlapping = FALSE, col="deepskyblue")
kpPlotRegions(kp, data=MS_translated_NLR_bed, avoid.overlapping = FALSE, col="red", data.panel=2)
# kpAddLabels(kp, "EG", cex=0.8, col="#888888")
# kpAddLabels(kp, "MS", data.panel=2, cex=0.8, col="#888888")

#Close the device
dev.off()


# EG
EG_genome_bed  <-
  read.delim("~/rstudio/bgimazia_musa/NLR_translated_genome/karyotype/Ensete_glaucum.genome.bed", header = T) %>%
  dplyr::rename(chr = seq.name) 
  # filter(!str_detect(chr, "Bscaffold")) %>%
  # mutate(chr = str_remove(chr, "B")) %>%
  # rbind(MS_genome_bed %>% tail(n=2))


EG_translated_NLR_bed <- read.delim("~/rstudio/bgimazia_musa/NLR_translated_genome/karyotype/Ensete_glaucum.tblastn.sorted.merged.bed", header = F)  %>%
  dplyr::rename(chr = V1,
                start  = V2,  
                end = V3)  %>%
  mutate(chr = str_remove(chr, "B"))

EG_predictedProteins_NLR_bed <- read.delim("~/rstudio/bgimazia_musa/NLR_predicted_proteins/karyotype/Ensete_glaucum.tblastn.sorted.merged.bed", header = F)  %>%
  dplyr::rename(chr = V1,
                start  = V2,  
                end = V3)  %>%
  mutate(chr = str_remove(chr, "B"))


## 
#Open the device
tiff("result_plots/EG_predicted_vs_translated_NLR_karyoplot.tiff", width = 800, height = 400)

pp <- getDefaultPlotParams(plot.type=2)
pp$data1outmargin <- 10
pp$data2outmargin <- 10
pp$topmargin <- 400
kp <- plotKaryotype(genome=EG_genome_bed, ideogram.plotter = NULL, plot.type=2, plot.params = pp, cex=1.5)
kpAddCytobandsAsLine(kp)
kpAddMainTitle(kp, "E. glaucum", cex=2)
kpPlotRegions(kp, data=EG_predictedProteins_NLR_bed, avoid.overlapping = FALSE, col="deepskyblue")
kpPlotRegions(kp, data=EG_translated_NLR_bed, avoid.overlapping = FALSE, col="red", data.panel=2)
# kpAddLabels(kp, "Predicted proteins", cex=0.8, col="#888888", side = "right")
# kpAddLabels(kp, "Translated genome", data.panel=2, cex=0.8, col="#888888", side = "right")

#Close the device
dev.off()


EG_translated_NLR_bed %>%
  group_by(chr) %>%
  summarise(count = n())


MA_predictedProteins_NLR_bed %>%
  group_by(chr) %>%
  summarise(MA = n()) %>%
  ungroup() %>%
  # select(MS_translated) %>%
  # sum()
  cbind(
    MB_predictedProteins_NLR_bed %>%
      group_by(chr) %>%
      summarise(MB = n()) %>%
      ungroup() %>%
      select(MB),
    
    MS_predictedProteins_NLR_bed %>%
      group_by(chr) %>%
      summarise(MS = n()) %>%
      ungroup() %>%
      select(MS),
    
    EG_predictedProteins_NLR_bed %>%
      group_by(chr) %>%
      summarise(EG = n()) %>%
      ungroup() %>%
      select(EG) %>%
      rbind(base::data.frame(EG= c("NA","NA")))
  ) %>%
  write.table("result_plots/NLR_loci_summary.txt", col.names = T, row.names = F, quote = F, sep = '\t')


# EV_mazia_translated_NLR_bed <-
  read.delim("~/rstudio/bgimazia_musa/NLR_translated_genome/karyotype/EV_mazia.tblastn.sorted.merged.bed", header = F)  %>%
  dplyr::rename(chr = V1,
                start  = V2,  
                end = V3)  %>%
    distinct(chr) %>% nrow()
  
  group_by(chr) %>%
  summarise(EV_mazia = n()) %>%
    ungroup() %>%
    select(EV_mazia) %>%
    sum()

# EV_mazia_predictedProteins_NLR_bed <- 
  read.delim("~/rstudio/bgimazia_musa/NLR_predicted_proteins/karyotype/EV_mazia.tblastn.sorted.merged.bed", header = F)  %>%
  dplyr::rename(chr = V1,
                start  = V2,  
                end = V3)  %>%
    distinct(chr) %>% nrow()
    group_by(chr) %>%
    summarise(EV_mazia = n()) %>%
    ungroup() %>%
    select(EV_mazia) %>%
    sum()

    
    
    # EV_mazia_translated_NLR_bed <-
    read.delim("~/rstudio/bgimazia_musa/NLR_translated_genome/karyotype/EV_bedadeti.tblastn.sorted.merged.bed", header = F)  %>%
      dplyr::rename(chr = V1,
                    start  = V2,  
                    end = V3)  %>%
      # distinct(chr) %>% nrow()
    
    group_by(chr) %>%
      summarise(EV_mazia = n()) %>%
      ungroup() %>%
      select(EV_mazia) %>%
      sum()
    
    # EV_mazia_predictedProteins_NLR_bed <- 
    read.delim("~/rstudio/bgimazia_musa/NLR_predicted_proteins/karyotype/EV_bedadeti.tblastn.sorted.merged.bed", header = F)  %>%
      dplyr::rename(chr = V1,
                    start  = V2,  
                    end = V3)  %>%
      # distinct(chr) %>% nrow()
    group_by(chr) %>%
      summarise(EV_mazia = n()) %>%
      ungroup() %>%
      select(EV_mazia) %>%
      sum()
    

library(tidyverse)
library(phylotools)

# import TEsorter results 
Species_List <- c( "EV_bedadeti", "EV_mazia", "Ensete_glaucum",
                  "Musa_acuminata", "Musa_balbisiana")

TEsorter_annotation <- c()


for (Each_Species in Species_List) {
  TEsorter_annotation <-
    read.delim(paste0("tesorter_out/",
                      # "tesorter_out/bedadeti.0.25_TE.rexdb-plant.dom.tsv"
                      Each_Species, "_TE.rexdb-plant.dom.tsv")) %>%
    mutate(genome = Each_Species) %>%
    rbind(TEsorter_annotation)

}

# reformat TEsorter result 
TEsorter_annotation.v1 <-
  TEsorter_annotation %>%
  mutate(TE = str_extract(X.id,"Class_I+.*[^|]"),
         TE = str_remove(TE,"\\:\\w+.*"),
         TE = str_remove(TE,"Class_\\w+\\/"),
         seq.name = str_remove(X.id,"Class_I+.*[^|]"),
         seq.name = str_remove(seq.name,"\\|"),
         seq.name = str_remove(seq.name,":\\d+-\\d+")) %>%
  select(length,evalue,coverge, genome, seq.name) %>%
  left_join(TEsorter_annotation_cls %>%
      mutate(seq.name = str_remove(X.TE,":\\d+-\\d+")) ) %>%
  select(seq.name,genome,length,coverge,evalue,Order,Superfamily,Clade) %>%
  mutate(TE = str_c(Order,Superfamily,Clade, sep = ':'))


# read predited proteins of EV mazia, Ev bedadeti and MA and MB, and execlude those with TE insertion in EV

Ensete_MUSA_seq <- c()

for (Species_List in Species_List){
  Ensete_MUSA_seq <- 
    read.fasta(paste0("gene_prediction_EDTA_masked/",Species_List,".fasta")) %>%
    mutate(genome = str_remove(Species_List,".maker.proteins.\\d+")) %>% 
    separate(seq.name, into = c("seq.name"), sep = ' ') %>%
    rbind(Ensete_MUSA_seq) 
}


  EV_prots = c("EV_mazia","EV_bedadeti")

  for (EV_prots_list in EV_prots) {
  Ensete_MUSA_seq %>%
    filter (str_detect(genome,"EV_\\w+")) %>%

    left_join(TEsorter_annotation.v1) %>% 
    mutate(length = str_replace_na(length,"NA")) %>%
    filter(length == 'NA') %>%
      filter(genome == EV_prots_list) %>%
      select(seq.name,seq.text) %>%
      dat2fasta(paste0("gene_prediction_EDTA_masked/",EV_prots_list,".TE_excluded.fasta"))
  }
      

# Summary 
    Ensete_MUSA_seq %>%
      left_join(TEsorter_annotation.v1) %>% 
      mutate(length = str_replace_na(length,"NA")) %>%
      filter(length == 'NA') %>%
      # filter(genome == Each_Species) %>%
      select(seq.name, genome) %>%
      mutate(seq.name = str_remove_all(seq.name,"-R.*")) %>%
      group_by(genome) %>%
      distinct() %>%
      summarise(count = n())

    

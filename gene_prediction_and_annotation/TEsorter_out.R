library(tidyverse)
library(ggtext)
library(phylotools)


# import TEsorter results 
Species_List <- c( "EV_bedadeti", "EV_mazia", "Ensete_glaucum",
                  "Musa_acuminata", "Musa_balbisiana")

TEsorter_annotation <- c()


TEsorter_annotation_cls <- c()

for (Each_Species in Species_List) {
  TEsorter_annotation <-
    read.delim(paste0("tesorter_out/",
                      # "gene_prediction/tesorter_out/bedadeti.0.25_TE.rexdb-plant.dom.tsv"
                      Each_Species, "_TE.rexdb-plant.dom.tsv")) %>%
    mutate(genome = Each_Species) %>%
    rbind(TEsorter_annotation)


  TEsorter_annotation_cls <-
    read.delim(paste0("tesorter_out/",
                      Each_Species, "_TE.rexdb-plant.cls.tsv")) %>%
    mutate(genome = Each_Species) %>%
    rbind(TEsorter_annotation_cls)


}



# read predited proteins of EV mazia, Ev bedadeti and MA and MB. 

Ensete_MUSA_seq <- c()

for (Species_List in Species_List){
  Ensete_MUSA_seq <- 
    read.fasta(paste0("gene_prediction_EDTA_masked/",Species_List,".fasta")) %>%
    mutate(genome = str_remove(Species_List,".maker.proteins.\\d+")) %>% 
    separate(seq.name, into = c("seq.name"), sep = ' ') %>%
    rbind(Ensete_MUSA_seq) 
}


predicted_prot_TEsorter_joined.AED <-
  Ensete_MUSA_seq %>% 
  mutate(seq.len = nchar(seq.text)) %>% 
  select(-seq.text) %>% 
           
  # joing find TE anntation
  left_join(
    TEsorter_annotation.v1) %>%
    mutate(coverage = round(length/seq.len,2),
           TE = str_replace_na(TE, "No_TE_found")) %>%
  filter(TE != "No_TE_found") 



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
      # rbind(read.delim("gene_prediction/tesorter_out/mazia_after_ncbi_submission_TE.rexd-pant.cls.tsv") %>%
      #         mutate(genome = "mazia")) %>%
      mutate(seq.name = str_remove(X.TE,":\\d+-\\d+")) ) %>%
  select(seq.name,genome,length,coverge,evalue,Order,Superfamily,Clade) %>%
  mutate(TE = str_c(Order,Superfamily,Clade, sep = ':'))
# 


TEsorter_annotation.v2 <-
  predicted_prot_TEsorter_joined.AED %>% 
  # filter(genome != "bedadeti.20", 
  #        genome != "bedadeti.30",
  #        genome != "bedadeti.35",
  #        genome != "bedadeti.40",
  #        #genome != "bedadeti.40",
  #         genome != "mazia.20",
  #        genome != "mazia.30",
  #        genome != "mazia.35",
  #        genome != "mazia.40",
  #        #genome != "mazia.40"
  #        ) %>%
  mutate(TE = str_replace(TE,"mixture\\:mixture","mixture"),
         TE = str_remove(TE,"\\:unknown"),
         genome = str_replace_all(genome, c( #"Musa_schizocarpa" = "MS", 
                                       "Musa_balbisiana" = "MB\nproteins", 
                                       "Musa_acuminata" = "MA\nproteins", 
                                       "EV_mazia" = "EV (Mazia)\nproteins (AED<=0.25)",  
                                       "Ensete_glaucum" = 'EG', 
                                       "EV_bedadeti" = "EV (Bedadeti)\nproteins (AED<=0.25)")),
         coverage = case_when(coverage >1 ~ 1,
                         TRUE ~ coverage))  


TEsorter_annotation.v2$genome <- factor(TEsorter_annotation.v2$genome,
                                        levels = c("EV (Mazia)\nproteins (AED<=0.25)", "EV (Bedadeti)\nproteins (AED<=0.25)",
                                                   "EG","MA\nproteins", "MB\nproteins",  "MS"))


# rearrange TE families for plotting 

TEsorter_annotation.v2$TE <- factor(TEsorter_annotation.v2$TE, 
                                    levels = c(
                                      "DIRS:unknown",
                                      "TIR:Tc1_Mariner",
                                      "TIR:Sola1",
                                      "TIR:PiggyBac",
                                      "TIR:PIF_Harbinger",
                                      "TIR:MuDR_Mutator",
                                      "TIR:hAT",
                                      "pararetrovirus:unknown",
                                      "mixture",
                                      "LINE:unknown",
                                      "LTR:mixture",
                                      "LTR:Copia:Bianca",
                                      "LTR:Copia",
                                      "LTR:Copia:Bryco",
                                      "LTR:Copia:Alesia",
                                      "LTR:Copia:Gymco-I",
                                      "LTR:Copia:Gymco-III",
                                      "LTR:Copia:Gymco-II",
                                      "LTR:Copia:Gymco-IV",
                                      "LTR:Copia:Lyco",
                                      "LTR:Copia:Osser",
                                      "LTR:Copia:TAR",
                                      "LTR:Copia:mixture",
                                      "LTR:Copia:Tork",
                                      "LTR:Copia:Ale",
                                      "LTR:Copia:Angela",
                                      "LTR:Copia:SIRE",
                                      "LTR:Copia:Ikeros",
                                      "LTR:Copia:Ivana",
                                      "LTR:Gypsy:Phygy",
                                      "LTR:Gypsy:Selgy",
                                      "LTR:Gypsy:Athila",
                                      "LTR:Gypsy:TatI",
                                      "LTR:Gypsy:chromo-outgroup",
                                      "LTR:Gypsy:chromo-unclass",
                                      "LTR:Gypsy:Chlamyvir",
                                      "LTR:Gypsy:TatII",
                                      "LTR:Gypsy:non-chromo-outgroup",
                                      "LTR:Gypsy:TatIII",
                                      "LTR:Gypsy:Tcn1",
                                      "LTR:Gypsy:Ogre",
                                      "LTR:Gypsy:mixture",
                                      "LTR:Gypsy:Tekay",
                                      "LTR:Gypsy:CRM",
                                      "LTR:Gypsy:Galadriel",
                                      "LTR:Gypsy:Reina",
                                      "LTR:Gypsy:Retand"))


## plot
TEsorter_annotation.v2 %>% 
  filter(genome != "EG",genome != "MS") %>%
  # filter(genome != "EG", genome != "MS") %>%
  ggplot(aes(TE,coverage, colour = genome))+
  geom_point(size = 3,  alpha = 0.5)+
  geom_jitter(width = 0.1, height = 0.1)+
  facet_grid(~genome)+
  coord_flip()+
  # scale_y_continuous(breaks =c(0.2,0.8))+
  scale_colour_manual( values = c("black","blue","green","red","gray","orange","#808000"))+
  # xlim(1,1000) +
  labs(y="Length coverage of (%) TE inserted genes in the predicted proteins",
       x = "Transposable elements\nOrder:Superfamily:Clade")+
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 9,face = "bold"),
    axis.title.y = element_text(size = 9, face = "bold",vjust = -8),
    axis.text.x = element_text(size = 8,face ="bold"),
    axis.text.y = element_text(size = 8,face = "bold"),
    strip.text.x = element_text(size=9, face = "bold"),
    legend.text = element_markdown(margin = margin(r=10),size = 10,face = "bold"),
    # legend.title = element_markdown(size=12,face = "bold"),
    legend.title = element_blank(),
    legend.box.background = element_rect(colour="black"),
    legend.position = "none",
    
  )

ggsave("tesorter_out/TE_sorter.EV_AED25.MAB.tiff", width=8, height=6)

# proportion of copy/gypsy 

TEsorter_annotation.v2 %>% 
  filter(genome != "EG",genome != "MS") %>% 
  filter(str_detect(Superfamily,"Gypsy") | str_detect(Superfamily,"Copia") ) %>% 
  group_by(genome,seq.name,Superfamily) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  select(-count) %>%
  group_by(genome,Superfamily) %>%
  summarise(count = n())


TEsorter_annotation.v2_TE_family_count <- 
TEsorter_annotation.v2 %>% 
  filter(genome != "EG",genome != "MS") %>% 
  filter(!str_detect(genome,"^M")) %>% 
  filter(Order !="mixture") %>%
  group_by(genome,TE) %>%
  summarise(count = n()) %>%
  mutate(TE = str_remove(TE,"LTR:"),
         TE = str_remove(TE,":unknown")) %>%
  # head()
  filter(TE !="mixture") 
  # head()
  
 TEsorter_annotation.v2_TE_family_count.1 <- 
TEsorter_annotation.v2_TE_family_count %>%
  mutate(type = str_extract(TE,"\\w+:"),
         TE_fam = str_remove(TE,"\\w+:"),
         type = str_replace_na(type, "NA")) %>%
  filter(type == "NA") %>%
  mutate(type = case_when(type == "NA" ~ TE)) %>%
  rbind(
    TEsorter_annotation.v2_TE_family_count %>%
      mutate(type = str_extract(TE,"\\w+:"),
             TE_fam = str_remove(TE,"\\w+:"),
             type = str_replace_na(type, "NA")) %>%
      filter(type != "NA")
  ) %>% 
  mutate(type = str_remove(type, ":")) 

  

# Plot figure showing number of predicted proteins with TE insertions for each TE families 


TEsorter_annotation.v2_TE_family_count.1 <- 
  TEsorter_annotation.v2_TE_family_count %>%
  mutate(type = str_extract(TE,"\\w+:"),
         TE_fam = str_remove(TE,"\\w+:"),
         type = str_replace_na(type, "NA")) %>%
  filter(type == "NA") %>%
  mutate(type = case_when(type == "NA" ~ TE)) %>%
  rbind(
    TEsorter_annotation.v2_TE_family_count %>%
      mutate(type = str_extract(TE,"\\w+:"),
             TE_fam = str_remove(TE,"\\w+:"),
             type = str_replace_na(type, "NA")) %>%
      filter(type != "NA")
  ) %>% 
  mutate(type = str_remove(type, ":")) 

TEsorter_annotation.v2_TE_family_count.1$TE_fam <- factor(TEsorter_annotation.v2_TE_family_count.1$TE_fam, 
                                                          levels = c(
                                                            "DIRS",
                                                            "Tc1_Mariner",
                                                            "Sola1",
                                                            "PiggyBac",
                                                            "PIF_Harbinger",
                                                            "MuDR_Mutator",
                                                            "hAT",
                                                            "pararetrovirus",
                                                            "LINE",
                                                            "mixture",
                                                            "Bianca",
                                                            "Copia",
                                                            "Bryco",
                                                            "Alesia",
                                                            "Gymco-I",
                                                            "Gymco-III",
                                                            "Gymco-II",
                                                            "Gymco-IV",
                                                            "Lyco",
                                                            "Osser",
                                                            "TAR",
                                                            # "mixture",
                                                            "Tork",
                                                            "Ale",
                                                            "Angela",
                                                            "SIRE",
                                                            "Ikeros",
                                                            "Ivana",
                                                            "Phygy",
                                                            "Selgy",
                                                            "Athila",
                                                            "TatI",
                                                            "chromo-outgroup",
                                                            "chromo-unclass",
                                                            "Chlamyvir",
                                                            "TatII",
                                                            "non-chromo-outgroup",
                                                            "TatIII",
                                                            "Tcn1",
                                                            "Ogre",
                                                            # "mixture",
                                                            "Tekay",
                                                            "CRM",
                                                            "Galadriel",
                                                            "Reina",
                                                            "Retand"))


TEsorter_annotation.v2_TE_family_count.1 %>% 
  ggplot(aes(TE_fam,count, colour = type))+
  geom_point(size = 4,  alpha = 5)+
  # geom_jitter(width = 0.1, height = 0.2)+
  facet_grid(~genome)+
  coord_flip()+
  scale_y_continuous(breaks =c(50,2500,5000,6746,8651))+
  # scale_colour_manual( values = c("black","blue","green","red","gray","orange","#808000"))+
  # xlim(1,1000) +
  labs(y="Number of TE enconding proteins families inserted in the predicted proteins ",
       x = "Transposable elements\nSuperfamily:Clade")+
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 9,face = "bold"),
    axis.title.y = element_text(size = 9, face = "bold",vjust = -9),
    axis.text.x = element_text(size = 8,face ="bold"),
    axis.text.y = element_text(size = 9,face = "bold"),
    strip.text.x = element_text(size=9, face = "bold"),
    legend.text = element_markdown(margin = margin(r=10),size = 10,face = "bold"),
    # legend.title = element_markdown(size=12,face = "bold"),
    legend.title = element_blank(),
    legend.box.background = element_rect(colour="black"),
    legend.position = c(0.85,0.8),
    legend.key.size = unit(0.1,"cm"),
    # legend.key.height = unit(0.1,"cm") 
    
  )
ggsave("tesorter_out/TE_sorter.EV_AED25.MAB.TE_families_count1.tiff", width=7, height=6)

# Frequeny plots 

# TE_removed_predicted_aa_length_frequency <- c()

## generate a frequency table 
make_freq_table <- function( lcl_list )
{
  ## This function will create a frequency table for 
  ## the one variable sent to it where that
  ## table gives the items, the frequency, the relative 
  ## frequeny, the cumulative frequency, the relative
  ## cumulative frequency, and the number of degrees to 
  ## allocate in a pie chart.
  ##
  ## The actual result of this function is a data frame 
  ## holding that table.
  lcl_freq <- table( lcl_list )
  lcl_size <- length( lcl_list )
  lcl_df <- data.frame( lcl_freq )
  names( lcl_df ) <- c("Items","Freq")
  # lcl_values <- as.numeric( lcl_freq )
  # lcl_df$rel_freq <- lcl_values /  lcl_size
  # lcl_df$cumul_freq <- cumsum( lcl_values )
  # lcl_df$rel_cumul_freq <- cumsum( lcl_values ) / lcl_size
  # lcl_df$pie <- round( 360*lcl_df$rel_freq, 1 )
  lcl_df
}


TE_frequency_per_gene <-  
    make_freq_table(TEsorter_annotation.v2 %>%
                      group_by(genome, seq.name) %>%
                      summarise(TE_count = n()) %>%
                      select(TE_count))

colnames(TE_frequency_per_gene) = c("Items", "Freq","pred_prot")

TE_frequency_per_gene %>% 
  filter(pred_prot !=0) %>% 
  ggplot(aes(Freq,pred_prot,color=Items))+
  geom_point(size=4)+
  scale_colour_manual( values = c("black","blue","green","red","gray","orange","#808000"))+
  scale_y_continuous(breaks = c(1,500,2500,4500, 6500, 9000)) +
  # scale_x_continuous(breaks =c(0.00, 0.25, 0.50, 0.75, 1.00))+
  labs(x= "Frequency of different TE classes\n inserted per single protein coding gene",
       y = " # of TE genes")+
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 16,face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12,face ="bold"),
    axis.text.y = element_text(size = 12,face = "bold"),
    strip.text.x = element_text(size=12, face = "bold"),
    legend.text = element_markdown(margin = margin(r=10),size = 12,face = "bold"),
    # legend.title = element_markdown(size=12,face = "bold"),
    legend.title = element_blank(),
    # legend.box.background = element_rect(colour="black"),
    legend.position = "top",
    
  )

ggsave("tesorter_out/TE_sorter.EV_AED40.MAB.frequency.tiff", width=5, height=5)


### summary tables 

TEsorter_annotation.v2 %>% 
  group_by(genome, seq.name) %>%
  summarise(TE_count = n()) %>%
  # head()
  ungroup() %>%
  group_by(genome) %>%
  summarise(count = n(), max = max(TE_count)) %>%
  write.table("tesorter_out/count_TE_inserted_gene.AED40.txt", col.names = T, row.names = F,
              quote = F, sep = '\t')
  

TEsorter_annotation.v2 %>% 
  group_by(genome, seq.name) %>%
  summarise(TE_count = n()) %>%
  ungroup() %>% 
  filter(TE_count >=1) %>% 
  select(-TE_count) %>%
  group_by(genome) %>%
  summarise(count=n())
  head()
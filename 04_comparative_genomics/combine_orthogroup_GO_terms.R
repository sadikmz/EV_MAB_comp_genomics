
#take EV_MAB.OG_GeneCount from OrthoFinder2 resulst

# load("scripts/OrthoFinder.EV.AED40.TE_removed.170.MABSIT.RData")

library(tidyverse)

EV_AED25_TE.OG_GeneCount %>% filter(Orthogroup == "OG0000332")

## combine GO terms of all speceis into one file 

Species_List <- c("mazia.AED25","bedadeti.AED25","Musa_acuminata",
                  "Musa_balbisiana")

# dir.create("../../../bgimazia_musa/reannotation_analysis/mazia_reannotation/orthofinder/go_terms/")

path_to_go_interpro = "/Users/u1866313/rstudio/bgimazia_musa/reannotation_analysis/mazia_reannotation/interproscan_out/"


GO_terms_EV_AED25_MAB <- c()

for (Each_Species in Species_List) {
  
  
  GO_terms_EV_AED25_MAB <- 
    read.delim(paste0(path_to_go_interpro, Each_Species, ".go_terms.tsv"), header = F) %>%
    rename(genes = V1,
           GO_ID = V2) %>%
    mutate(genome = Each_Species) %>%
    rbind(GO_terms_EV_AED25_MAB)
}


GO_terms_EV_AED25_MAB %>% 
  filter(str_detect(GO_ID,"^GO")) %>%
  tail()

# split multiple GO terms of a gene to column data   
GO_terms_EV_AED25_MAB_v1 <-
  GO_terms_EV_AED25_MAB %>% 
  filter(GO_ID != "", GO_ID != "") %>%
  filter(str_detect(GO_ID,"^GO")) %>%
  distinct() %>% 
  mutate(genome = str_replace_all(genome, c("Ensete_ventricosum_bedadeti"="bedadeti.AED25",
                                            "Ensete_ventricosum_mazia"="mazia.AED25")),
    sep_GO_ID = str_count(GO_ID,"\\|")) %>% 
  # filter(sep_GO_ID > 7) %>% head()
  mutate(
    comb = str_c(genes,genome,sep_GO_ID,GO_ID, sep = "|")
  ) %>% 
  select(comb) %>% 
  separate(comb, c("genes","genomes", "OG_count", "GO_ID_set1","GO_ID_set2","GO_ID_set3",
                   "GO_ID_set4","GO_ID_set5","GO_ID_set6",
                   "GO_ID_set7","GO_ID_set8"), sep = "\\|") 



#combine into one file 

GO_terms_EV_AED25_MAB_v2 <-
  GO_terms_EV_AED25_MAB_v1 %>% 
  distinct() %>%
  select(genes,genomes,OG_count,GO_ID_set1) %>%
  bind_rows(
    GO_terms_EV_AED25_MAB_v1 %>%
      select(genes,genomes,OG_count,GO_ID_set2) %>% 
      rename(GO_ID_set1 = GO_ID_set2),
    GO_terms_EV_AED25_MAB_v1 %>%
      select(genes,genomes,OG_count,GO_ID_set3) %>% 
      rename(GO_ID_set1 = GO_ID_set3),
    GO_terms_EV_AED25_MAB_v1 %>%
      select(genes,genomes,OG_count,GO_ID_set4) %>% 
      rename(GO_ID_set1 = GO_ID_set4),
    GO_terms_EV_AED25_MAB_v1 %>%
      select(genes,genomes,OG_count,GO_ID_set5) %>% 
      rename(GO_ID_set1 = GO_ID_set5),
    GO_terms_EV_AED25_MAB_v1 %>%
      select(genes,genomes,OG_count,GO_ID_set6) %>% 
      rename(GO_ID_set1 = GO_ID_set6),
    GO_terms_EV_AED25_MAB_v1 %>%
      select(genes,genomes,OG_count,GO_ID_set7) %>% 
      rename(GO_ID_set1 = GO_ID_set7),
    GO_terms_EV_AED25_MAB_v1 %>%
      select(genes,genomes,OG_count,GO_ID_set8) %>% 
      rename(GO_ID_set1 = GO_ID_set8),
  ) %>%
  rename(GO_terms  = GO_ID_set1) %>%
  mutate( GO_terms = str_replace_na(GO_terms,"No_GO_terms")) %>% 
  select(-OG_count) %>%
  filter( GO_terms != 'No_GO_terms') 


GO_terms_EV_AED25_MAB_v2 %>% head()

# merge all genes with with out with out GO terms 

# bedadeti and mazia predicted proteins AED=0-0.40 and aa lenth >= 70 softlinked to working directory from blasp direcoty 
library(phylotools)

# EV_EG_MUSA = c("musa_ac","musa_ba","mazia","bedadeti")

# EV_EG_MABS_prots_ID_len <- c()
# 
# for (Each_MUSA_EG in EV_EG_MUSA) {
#   EV_EG_MABS_prots_ID_len <- 
# read.fasta(paste0("../musa_predicted_protiens/", Each_MUSA_EG, ".fasta")) %>%
#   mutate(prot_len = nchar(seq.text),
#          genome = Each_MUSA_EG) %>%
#     select(seq.name, prot_len, genome) %>%
#     rbind(
#       EV_EG_MABS_prots_ID_len
#     )
# }



Ensete_MUSA_seq_TE <-
  read.delim("~/rstudio/bgimazia_musa/reannotation_analysis/mazia_reannotation//Ensete_MUSA_seq.txt") %>% 
  filter(genome != "mazia", genome != "bedadeti", genome != "Ensete_glaucum", genome != "Musa_schizocarpa") %>% 
  # rbind(
  #   read.fasta("~/rstudio/bgimazia_musa/reannotation_analysis/mazia_reannotation/gene_prediction/gene_models/bedadeti.AED.0-0.40.non_TE.fasta") %>%
  #     mutate(genome = "bedadeti"),
  #   read.fasta("~/rstudio/bgimazia_musa/reannotation_analysis/mazia_reannotation/gene_prediction/gene_models/mazia.AED.0-0.40.non_TE.fasta") %>%
  #     mutate(genome = "mazia") 
  # ) 
  rbind(read.fasta("~/rstudio/bgimazia_musa/reannotation_analysis/mazia_reannotation/TE_inserted_prots/bedadeti.longest_isoforms.fasta") %>%
          mutate(genome = "bedadeti"),
        read.fasta("~/rstudio/bgimazia_musa/reannotation_analysis/mazia_reannotation/TE_inserted_prots/mazia.longest_isoforms.fasta") %>%
          mutate(genome = "mazia")) #%>%colnames()
# select(-gene_id,-length))

# Ensete_MUSA_seq %>% distinct(genome)

GO_terms_EV_AED25_TE <-
  Ensete_MUSA_seq_TE %>%
  select(-seq.text) %>% 
  # EV_EG_MABS_prots_ID_len %>% 
  # filter(genome == "musa_it") %>% 
  # mutate(seq.name.mit = str_extract(seq.name,"\\w+[^\\s]"),
  #        seq.name = case_when(str_detect(seq.name.mit,"^\\Mi_") ~ seq.name.mit,
  #                             TRUE ~ seq.name)) %>% 
  # dplyr::select(-seq.name.mit) %>% 
  left_join(
    GO_terms_EV_AED25_MAB_v2 %>%
      dplyr::rename(seq.name = genes) ) %>% 
  mutate(GO_terms = str_replace_na(GO_terms, 'NA'),
         genomes = str_replace_na(genomes,"NA" ),
         genomes = str_replace_all(genomes,c("mazia.AED25"="Ensete_ventricosum_mazia",
                                             "bedadeti.AED25"="Ensete_ventricosum_bedadeti")),
         genomes = str_replace_na(genomes,"NA" )) %>%
  # filter(GO_terms != "NA") %>%
  # filter(genomes == "NA") %>%
  mutate( genomes = case_when(genome == "bedadeti" & genomes == "NA" ~"Ensete_ventricosum_bedadeti",
                              genome == "mazia" & genomes == 'NA' ~ "Ensete_ventricosum_mazia",
                              # genome == "ensete_gl" & genomes == 'NA' ~ "Ensete_glaucum",
                              # genome == "musa_sc" & genomes == 'NA' ~ "Musa_schizocarpa",
                              genome == "musa_ba" & genomes == 'NA' ~ "Musa_balbisiana",
                              genome == "musa_ac" & genomes == 'NA' ~  "Musa_acuminata",
                              genome == "NA" & genomes == 'Musa_acuminata' ~  "Musa_acuminata",
                              genome == "Musa_acuminata" & genomes == 'NA' ~  "Musa_acuminata",
                              genome == "Musa_balbisiana" & genomes == 'NA' ~  "Musa_balbisiana",
                              genome == "NA" & genomes == 'Musa_balbisiana' ~  "Musa_balbisiana",
                              TRUE ~ genomes)) 


# GO terms with associated description 

# BiocManager::install("GO.db")
library(GO.db)

goterms = unlist(Term(GOTERM))
goterms_table <-
  goterms %>% as.data.frame() %>%
  rownames_to_column("GO_terms")
# 
colnames(goterms_table) = c("GO_terms", "description")

library(tidyverse)

goterms_table %>% head()
GO_terms_EV_AED25_TE %>% filter(GO_terms != "NA") %>% head()


## attached GO temrs with gene id
GO_terms_EV_AED25_TE_MAB_desc <-
  GO_terms_EV_AED25_TE %>%
  distinct() %>%
  left_join(
    goterms_table
  ) %>%
  mutate(description = str_replace_na(description,"NA"))

# number of genes that have go terms 
# GO_terms_EV_AED25_MAB_desc %>%
#   # filter(genomes == "Musa_balbisiana") %>%
#   mutate(genes = str_remove_all(seq.name,"-R\\w+")) %>% 
#   filter(str_detect(GO_terms,"GO:\\d+")) %>% 
#   dplyr::select(genomes, genes) %>%
#   distinct() %>% 
#   group_by(genomes) %>%
#   summarise(count = n()) %>%
#   write.table("go_terms/number_of_genes_with_GO_terms_EV_AED25_MAB", sep = '\t', col.names = T, quote = F,row.names = F)


# GO_terms_EV_AED25_MAB_v3 <-
# GO_terms_EV_AED25_MAB_v2 %>% 
#   distinct() %>%
#   group_by(genomes, GO_terms) %>% 
#   summarise( count = n()) %>% 
#   # filter( count >= 5000)
#   left_join(
#     goterms_table
#   ) %>% 
# mutate(genomes = str_replace_all(genomes,"Ensete_ventricosum_bedadeti", "EV (Bedadeti)"),
#        genomes = str_replace_all(genomes,"Ensete_ventricosum_mazia", "EV (Mazia)"),
#        genomes = str_replace_all(genomes, "EG"),
#        genomes = str_replace_all(genomes, "Musa_balbisiana", "MB"),
#        genomes = str_replace_all(genomes, "Musa_acuminata", "MA"),
#        genomes = str_replace_all(genomes, "MS")) 
# group_by(genome, GO_terms) %>%
# summarise(count=n()) %>% 
# filter(count >= 1) %>% head()
# ggplot(aes(count)) +
# geom_freqpoly()
# nrow()


## genrate GO_terms labels 
GO_terms_EV_AED25_TE_MAB_desc.v1 <-
  GO_terms_EV_AED25_TE_MAB_desc %>% 
    # separate(seq.name,into = "seq.name", sep = ' ') %>%
  distinct() %>%
  mutate(genomes_short = str_replace_all(genomes, c("Ensete_ventricosum_bedadeti"="EV (Bedadeti)",
                                                    "Ensete_ventricosum_mazia"="EV (Mazia)", 
                                                    # "Ensete_glaucum"= "EG",
                                                    "Musa_balbisiana" = "MB", 
                                                    "Musa_acuminata" ="MA"
                                                    # "Musa_schizocarpa" = "MS"
  )))

GO_terms_EV_AED25_TE_MAB_desc.v1 %>% tail(n=10)

# GO_terms_EV_AED25_MAB_desc.v1 <-   
# GO_terms_EV_AED25_MAB_desc.v1$genomes_short = factor(GO_terms_EV_AED25_MAB_desc.v1$genomes_short, levels = c("EV (Bedadeti)","EV (Mazia)",
# "EG","MB","MA","MS")) 


## combine GO_term with Orthogroups parse gene list: all_assigned_OGsgenes

# assigned_unassigned_OGs_gene_EV_EG_MABS <- 
#   all_assigned_OGsgenes_EV_EG_MABS %>% 
#   dplyr::rename(seq.name = pred_genes,
#                 genomes = species) %>% 
#   rbind(Orthogroups_UnassignedGenes_EV_EG_MBS %>% 
#           dplyr::rename(genomes = genome)) 
# 
# GO_terms_EV_AED25_MAB_OGs_genes <-
#   GO_terms_EV_AED25_MAB_desc.v1 %>% 
#   left_join(
#     assigned_unassigned_OGs_genes_EV_EG_MAB %>%
#       mutate(genome = str_replace_all(genome,c("musa_ac"="Musa_acuminata",
#                                                "musa_ba" = "Musa_balbisiana")))
#   ) %>%
#   mutate(Orthogroup = str_replace_na(Orthogroup,"")) 

GO_terms_EV_AED25_OGs_genes
GO_terms_EV_AED25_MAB_OGs_genes <-
  assigned_unassigned_OGs_genes_TE_EV %>%
  # mutate(genome = str_replace_all(genome,c("musa_ac"="Musa_acuminata",
  #                                          "musa_ba" = "Musa_balbisiana"))) %>%
  left_join(GO_terms_EV_AED25_TE_MAB_desc.v1 %>%
              filter(!str_detect(genome,"^Musa")) #%>%
              # mutate(genome = str_replace_all(genome,c("mazia"="Ensete_ventricosum_mazia",
              #                                          "bedadeti" = "Ensete_ventricosum_bedadeti")))
              ) %>%
  mutate(GO_terms  = str_replace_na(GO_terms,"NA")) %>%
    select(-genomes,-genomes_short)  %>%
  filter(GO_terms !="NA") 
  

  
#   
#   filter(genome == "musa_ac", str_detect(Orthogroup, "^OG")) %>%
#   filter(str_detect(seq.name,"Macma4_Mt_g00010.1"))
# #   distinct(genome)
# 
# GO_terms_EV_AED25_MAB_OGs_genes %>% filter(genomes=="NA")
GO_terms_EV_AED25_MAB_OGs_genes  %>% 
  filter(genome == "Musa_balbisiana") %>%
  # filter(GO_terms !="NA") %>%
  arrange(desc(Orthogroup)) %>%
  nrow()


GO_terms_EV_AED25_TE_MAB_desc.v1 %>%
  filter(genome=="bedadeti") %>%
  select(seq.name) %>%
  distinct() %>% 
  nrow()


GO_terms_EV_AED25_TE_MAB_desc.v1 %>%
  filter(genome=="mazia") %>%
  select(seq.name) %>%
  distinct() %>% 
  # head(n=20)
  nrow()



assigned_unassigned_OGs_genes_EV_EG_MAB %>%
  filter(genome=="Ensete_ventricosum_bedadeti") %>%
  distinct() %>%
  head()

assigned_unassigned_OGs_genes_EV_EG_MAB %>%
  filter(genome=="Ensete_ventricosum_mazia") %>%
  distinct() %>%
  nrow()



save.image(file = "../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED25_MAB/OrthoFinder.EV.AED25.TE.RData")
# load(file = "../go_terms/GO_terms_EV170_EG_MAB.RData")



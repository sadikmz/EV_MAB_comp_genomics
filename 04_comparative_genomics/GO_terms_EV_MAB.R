# This script parse GO_terms annotation generated for predicted protein coding genes and link it with transcript ID of each gene. 

# List of package required
packages = c("tidyverse", "phylotools")

# Install packages missing packages.
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}


library(tidyverse)
library(phylotools)


## combine GO terms of all speceis into one file 

Species_List <- c("mazia.AED25","bedadeti.AED25","Musa_acuminata",
                  "Musa_balbisiana","Ensete_glaucum")

dir.create("absoute_path_to_orthofinder/go_terms/")

path_to_go_interpro = "absolute_path_to_interproscan_out/"


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

# a combined predicted proteins of mazia and bedadeti, Ensete glaucum, and Musa species.

library(phylotools)

Prot_List = c("Ensete_glaucum",  "Musa_acuminata","Musa_balbisiana","mazia.AED25","bedadeti.AED25")

Ensete_MUSA_seq <- c()

for (K in Prot_List) {
  Ensete_MUSA_seq <-
read.fasta(paste0("data_pred_prots/", K, ".fasta")) %>%
  mutate(prot_len = nchar(seq.text),
         genome = Prot_List) %>%
    rbind(Ensete_MUSA_seq)
}


# Ensete_MUSA_seq %>% distinct(genome)

GO_terms_EV_AED25_MAB_v3 <-
  Ensete_MUSA_seq %>%
  select(-seq.text) %>% 
  left_join(
    GO_terms_EV_AED25_MAB_v2 %>%
      dplyr::rename(seq.name = genes) ) %>% 
  mutate(GO_terms = str_replace_na(GO_terms, 'NA'),
         genome = str_replace_na(genome,"NA" ),
         genome = str_replace_all(genome,c("mazia.AED25"="Ensete_ventricosum_mazia",
                                             "bedadeti.AED25"="Ensete_ventricosum_bedadeti")),
         genome = str_replace_na(genome,"NA" )) 


# Add updated descriptions of GO-terms to IDs from GO.db that contain annotation maps describing the entire Gene Ontology.

BiocManager::install("GO.db")
library(GO.db)

goterms = unlist(Term(GOTERM))
goterms_table <-
  goterms %>% as.data.frame() %>%
  rownames_to_column("GO_terms")
# 
colnames(goterms_table) = c("GO_terms", "description")

library(tidyverse)

goterms_table %>% head()
GO_terms_EV_AED25_MAB_v3 %>% filter(GO_terms != "NA") %>% head()


## attached GO temrs with gene id
GO_terms_EV_AED25_MAB_desc <-
  GO_terms_EV_AED25_MAB_v3 %>%
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
GO_terms_EV_AED25_MAB_desc.v1 <-
  GO_terms_EV_AED25_MAB_desc %>% 
    # separate(seq.name,into = "seq.name", sep = ' ') %>%
  distinct() %>%
  mutate(genomes_short = str_replace_all(genomes, c("Ensete_ventricosum_bedadeti"="EV (Bedadeti)",
                                                    "Ensete_ventricosum_mazia"="EV (Mazia)", 
                                                    # "Ensete_glaucum"= "EG",
                                                    "Musa_balbisiana" = "MB", 
                                                    "Musa_acuminata" ="MA"
                                                    # "Musa_schizocarpa" = "MS"
  )))


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


# GO_terms_EV_AED25_MAB_OGs_genes <-
#   assigned_unassigned_OGs_genes_EV_EG_MAB %>%
#   mutate(genome = str_replace_all(genome,c("musa_ac"="Musa_acuminata",
#                                            "musa_ba" = "Musa_balbisiana"))) %>%
#   left_join(GO_terms_EV_AED25_MAB_desc.v1 %>%
#               mutate(genome = str_replace_all(genome,c("mazia"="Ensete_ventricosum_mazia",
#                                                        "bedadeti" = "Ensete_ventricosum_bedadeti")))) %>%
#   mutate(GO_terms  = str_replace_na(GO_terms,"NA")) %>%
#     select(-genomes,-genomes_short)  %>%
#   filter(GO_terms !="NA") 
  

# save
save.image(file = "OrthoFinder.EV.AED25.MAB.RData")



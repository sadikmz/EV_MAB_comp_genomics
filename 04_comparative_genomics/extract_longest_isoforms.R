# The script extract the longest predicted protein isoforms for orthology analysis. 
# It assumes a maximum of two equal length isoforms and it need be modified in cases where more than three equal length isoforms of a single gene exist.

# Load library


# List of package required
packages = c("tidyverse", "phylotools", "janitor")

# Install packages missing packages.
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

library(tidyverse)
library(phylotools)
library(janitor)

# list.files("gene_prediction_EDTA_masked/", pattern = ".fasta")

# input_dir = "path_to_input_directory"

extract_longest_protein_isoform <-
  
function(input_dir,
         input_fasta,
         output_file,
         geneID_pattern){
  
  input_fasta <- file.path(input_dir, basename(input_fasta))
  
  
# bedadeti
input_prot <-
  phylotools::read.fasta(input_fasta)  %>% 
  mutate(len = nchar(seq.text),
         gene.name = str_extract(seq.name, geneID_pattern)
         # seq.name= str_remove(seq.name,"protein.*"),
  ) %>% 
  # filter(len >=50) %>%
  distinct()



# extract longest protein isoforms 
longest_isoform <-
  input_prot %>%
  group_by(gene.name) %>%
  summarise(len=max(len)) 

# identifify equal length isoforms

uncleaned_longest_isoform_fasta <-
  longest_isoform %>% 
  select(gene.name) %>% duplicated() %>%  
  as_tibble() %>% 
  bind_cols(longest_isoform) %>% 
  filter(value == "FALSE")  %>% # if there are equale lenght isoforms, 
  # check their aa length correspondes to T/F assignment in value parameter. 
  mutate(value = str_replace_na(value, "TRUE")) %>% 
  # join input fasta
  left_join(
    input_prot
  ) %>% 
  distinct() 


duplicated_seq.name <- 
  uncleaned_longest_isoform_fasta  %>%
  select(gene.name) %>%
  get_dupes() %>%
  distinct(gene.name) 


## equal length isoforms 

isoforms_count <- length(duplicated_seq.name$gene.name)

equal_length_isoforms <-
  uncleaned_longest_isoform_fasta %>% 
  left_join(
    duplicated_seq.name %>%
      mutate(remove = "T")
  ) %>% 
  mutate(remove = str_replace_na(remove, "F")) %>% 
  filter(remove == "T" ) %>% 
  bind_cols(
    data.frame( select= rep(c(1,2),isoforms_count))
  ) %>%
  filter(select == 1)  %>%
  select(-select, -remove, -len, -value)


# longest_isoform_fasta
longest_isoform_fasta <-
  uncleaned_longest_isoform_fasta %>% 
  left_join(
    duplicated_seq.name %>%
      mutate(remove = "T")
  ) %>% 
  mutate(remove = str_replace_na(remove, "F")) %>% 
  filter(remove == "F") %>% 
  # filter(seq.name =="EVEPO072292")
  distinct() %>%
  select(-remove,-value, -len) %>% 
  bind_rows (
    equal_length_isoforms
  )

longest_isoform_fasta %>% 
  separate(seq.name,into = "seq.name", sep = ' ') %>%
  select(seq.name,seq.text) %>% 
  dat2fasta(paste0(input_dir, output_file,".longest_isoforms.fasta"))
}

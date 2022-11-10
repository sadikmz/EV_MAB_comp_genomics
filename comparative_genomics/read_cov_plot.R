
library(tidyverse)

bed_cov_list=list.files(path="align_out", pattern="coverage.bed")

datalist <- ()


for (i in bed_cov_list){
  dat <- read.delim(paste0("align_out/",i, ".coverage.bed"),header = F) %>%
  rename(chr=V1,
         start = V2,
         end = V3,
         coverage = V4)

  dat$genome <- i
  datalist[[i]] <- dat # add it to your list
}

EV_MAB_mapping_coverage = do.call(rbind, datalist)

save.image(file = "mazia_EG_MABS_alignment_coverage_dedup.RData")

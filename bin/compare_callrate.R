#!/usr/bin/env Rscript
# For pairs of genetic duplicate genotype samples, find those with non-optimal
# call rate and prepare file suitable to remove them with plink

library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
file_imiss  <- args[1]
file_dups <- args[2]
file_out <- args[3]

callrate <- fread(file_imiss)
dups.same_IID <- fread(file_dups)

callrate %>% 
  filter(IID %in% c(dups.same_IID$V2, dups.same_IID$V4)) %>% 
  mutate(IID_clean = str_remove(IID, "_.+$")) %>% 
  group_by(IID_clean) %>% mutate(best_callrate = F_MISS == min(F_MISS)) %>%
  group_by(IID_clean, F_MISS) %>% arrange(IID) %>% 
  mutate(remove = !best_callrate | row_number() != 1) %>%
  filter(remove) %>% ungroup() %>%
  select(FID, IID) %>%
  fwrite(file_out, col.names = FALSE, sep = "\t")
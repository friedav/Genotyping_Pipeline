#!/usr/bin/env Rscript
# Update variant IDs in BIM file by converting Illumina IDs to rsIDs based on
# mapping file provided by Illumina for each of their arrays;
# if this introduces duplicated variants, write table with variants to remove

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
file.bim <- args[1]
file.idmap <- args[2]

# Illumina mapping table:
# - in case multiple rs IDs are mapped to one Illumina ID (e.g. "rs1,rs3,rs4" in
#   the RsID column), consider only the first one (e.g. "rs1")); 
# - remove rows without mapped rs ID, i.e. RsID == "."
idmap <- fread(file.idmap, data.table = FALSE) %>% 
  mutate(RsID = str_remove(RsID, ",.+$")) %>% 
  filter(RsID != ".")

# target bim file
bim <- fread(file.bim, col.names = c("chr", "id", "cm", "pos", "a1", "a2")) 

# update non-rs IDs (if present in mapping table) 
bim$id <- with(bim, ifelse(grepl("^rs[0-9]+$", id) | !id %in% idmap$Name, 
                           id, 
                           idmap$RsID[match(id, idmap$Name)]))

# identify and flag duplicated variants
bim <- bim %>% group_by(id) %>% 
  mutate(idx = 1:n(),
         id = ifelse(idx == 1, id, paste0(id, "_dupl", idx-1))) %>% 
  ungroup() %>% select(-idx)

# write duplicate exclusion file to use with Plink
exclude <- grep("dupl", bim$id, value = TRUE)
write(exclude, str_replace(file.bim, ".bim$", ".duplicates.exclude"))

fwrite(bim, str_replace(file.bim, ".bim$", ".updatedIDs_dupsFlagged.bim"),
       col.names = FALSE, sep = "\t")

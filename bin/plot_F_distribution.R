#!/usr/bin/env Rscript
# Plot distribution of F estimates, either for autosomal heterozygosity or
# x-chromosomal heterozygosity colored by sex;
# in case of autosomal heterozygosity, F estimates are standardized to mean 0
# and sd 1
#
# Arguments:
# 1. Path to Plink output with F estimates (*.het or *.sexcheck)
# 2. Optional: comma-sep string of thresholds to add to plot (in case of *.het,
#    the values are interpreted as standard deviations from mean)

library(ggplot2)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
file_fhet  <- args[1]
thres <- args[2]

fhet <- fread(file_fhet)

# x-chr het if PEDSEX and SNPSEX columns are present
if ("PEDSEX" %in% colnames(fhet) & "SNPSEX" %in% colnames(fhet)) {
  plt <- ggplot(fhet, aes(F, fill = as.factor(PEDSEX))) + 
      geom_histogram(bins = 50)
  
} else {
  plt <- ggplot(fhet, aes(scale(F))) + geom_histogram(bins = 50)
}

# add threshold to plot if provided
if (!is.na(thres)) {
  thres <- as.numeric(unlist(strsplit(thres, ",")))
  plt <- plt + geom_vline(xintercept = thres)
}

ggsave(paste0(file_fhet, ".Fdistrib.png"), plot = plt,
       width = 15, height = 12, dpi = 600)

#!/usr/bin/env Rscript
# Plot distribution of F estimates by sex

library(ggplot2)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
file_sexcheck  <- args[1]

sexcheck <- fread(file_sexcheck)

plt <- ggplot(sexcheck, aes(F, fill = as.factor(PEDSEX))) + 
    geom_histogram(bins = 50)
ggsave(paste0(file_sexcheck, ".Fdistrib.png"), plot = plt,
       width = 15, height = 12, dpi = 600)


#!/usr/bin/env Rscript
#
# Identify genetic outliers: distance from the mean of >x*SD in the given number
# of ancestry components from principal component analysis (PCA)

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
thres <- as.integer(args[1])
file.pca <- args[2]
outpref <- args[3]

pca <- fread(file.pca) %>% mutate(outlier = NA_character_) %>% 
  mutate(across(matches("PC[0-9]+"), scale))

# start from highest dimension to report lowest dimension with outlier status
for (comp in grep("PC[0-9]+", colnames(pca), value = TRUE) %>% rev()) {
  pca <- mutate(pca, outlier = ifelse(abs(pca[[comp]]) > thres, comp, outlier))
}

# avoid ggplot error: "Must request at least one colour from a hue palette."
if (all(is.na(pca$outlier))) pca$outlier <- "None"

# plot PC1 and PC2
plot <- ggplot(pca, aes(PC1, PC2, col = outlier)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = c(-thres, thres), linetype = "dotted", color = "red") +
  geom_hline(yintercept = c(-thres, thres), linetype = "dotted", color = "red") +
  geom_point(alpha = 0.5) +
  theme_minimal()
ggsave(paste0(outpref, "_plot.pdf"), plot = plot)


# write list of pca outliers if any
if (any(!is.na(pca$outlier)) & !all(pca$outlier == "None")) {
  file.out <- paste0(outpref, ".remove")
  pca %>% filter(!is.na(outlier)) %>% select(FID, IID, outlier) %>% 
    fwrite(file.out, sep = "\t")
}

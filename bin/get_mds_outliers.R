#!/usr/bin/env Rscript
#
# Identify genetic outliers: distance from the mean of >x*SD in the given number
# of multidimensional scaling (MDS) ancestry components

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
thres <- as.integer(args[1])
file.mds <- args[2]
outpref <- args[3]

mds <- fread(file.mds) %>% mutate(outlier = NA_character_) %>% 
  mutate(across(matches("C[0-9]+"), scale))

# start from highest dimension to report lowest dimension with outlier status
for (comp in grep("C[0-9]+", colnames(mds), value = TRUE) %>% rev()) {
  mds <- mutate(mds, outlier = ifelse(abs(mds[[comp]]) > thres, comp, outlier))
}

# avoid ggplot error: "Must request at least one colour from a hue palette."
if (all(is.na(mds$outlier))) mds$outlier <- "None"

# plot C1 and C2
plot <- ggplot(mds, aes(C1, C2, col = outlier)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = c(-thres, thres), linetype = "dotted", color = "red") +
  geom_hline(yintercept = c(-thres, thres), linetype = "dotted", color = "red") +
  geom_point(alpha = 0.5) +
  theme_minimal()
ggsave(paste0(outpref, "_plot.pdf"), plot = plot)


# write list of mds outliers if any
if (any(!is.na(mds$outlier)) & !all(mds$outlier == "None")) {
  file.out <- paste0(outpref, ".remove")
  mds %>% filter(!is.na(outlier)) %>% select(FID, IID, outlier) %>% 
    fwrite(file.out, sep = "\t")
}

#!/usr/bin/env Rscript
#
# Swap lines in Plink FAM file, e.g. to resolve sample swaps
#
# If sample swaps in a Plink genotype bfile set have been identified, they can
# be resolved by swapping the respective lines in the FAM file, since the sample
# IDs are assigned to the actual genotype data in the BED file via their
# position. Note: By swapping full lines including sex info, this might also
# resolve sex mismatches detected during data quality control.
#
# Usage:
#     swap_samples.R \
#         <input fam file> \
#         <swapping file> \
#         <output fam file>
#
# Swapping file: Tab-separated text file with the old sample IID in the 1st and
#                the new sample IID in the 2nd column; this also allows for
#                triangle swaps; must have column names "current" and
#                "replacement"

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
file_fam_in  <- args[1]
file_swap    <- args[2]
file_fam_out <- args[3]

fam <- fread(file_fam_in, header = FALSE,
             col.names = c("FID", "IID", "PID", "MID", "sex", "pheno"))
swap <- fread(file_swap)

# lines containing sample info that will be replaced / will act as replacement
idx_current <- match(swap$current, fam$IID)
idx_replacement <- match(swap$replacement, fam$IID)

# check if old and new sample IIDs were actually found in fam file
if (any(is.na(c(idx_current, idx_replacement)))) stop("Error: Check swapping file.")

replacement <- fam[idx_replacement, ]
fam[idx_current, ] <- replacement

fwrite(fam, file = file_fam_out, sep = " ", col.names = FALSE)

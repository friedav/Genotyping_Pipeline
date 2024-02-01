#!/usr/bin/env bash
#
# Quality control and preprocessing of genotype data from PLINK file set: Part 2
#
# The preimputation genotype QC consists of two parts, since in case of the
# identification of sex mismatches and genetic duplicates in part 1 you might
# want to investigate and resolve some of these issues before continuing with
# part 2, instead of simply excluding all these samples.
#
# Prerequesits: Sex column in FAM file reflects phenotypic sex info if available
#
# Usage:        
# 02_run_genotype_QC_part2.sh \
#       <prefix used for bfile sets of processed data> \
#       <output directory> \
#       <Illumina variant ID to rs ID mapping file (array-specific)> \
#       <optional: remove file> \
#       <optional: swap file> 
#
# General thoughts: 
# 1. Is technical data quality okay (variants/individuals)? -> Call rate, Fhet
# 2. Is individual identity as expected? -> Sex mismatches, genetic duplicates
set -e

# mamba env create -f environment.yml -p ./env
conda activate ./env

outpref=$1
outdir=$2
out=$outdir/$outpref

illumina2rsid=$3

remove=$4
swaps=$5

# auxiliary data files / directories
king_1kg_ref=data/king_1kg_ref/KGref



################################################################################
#### Resolve sample swaps and other sample issues, if any
################################################################################

# sample swaps that can be resolved by swapping back IDs
if [ -z "$swaps" ]; then
    echo No sample swaps are resolved
    famresolved="${out}_QC4.fam"
else
    echo Resolving sample swaps: 
    cat $swaps
    famresolved=${out}_QC4.resolved_swaps.fam 
    Rscript bin/swap_samples.R ${out}_QC4.fam ${swaps} ${famresolved}
fi

# other sample identity issues identified via sex mismatches and genetic dups
if [ -z "$remove" ]; then
    cat ${out}_QC4_dups.notBestCallrate.remove | sort | uniq > ${out}_QC4.remove
else
    sed 's/$/\tSample identity (genetic duplicate and or sex mismatch)/' \
        ${remove} >> ${out}.sample_exclusions.tsv
    
    cat ${remove} \
        ${out}_QC4_dups.notBestCallrate.remove | sort | uniq > ${out}_QC4.remove
fi

plink --bfile ${out}_QC4 --fam ${famresolved} --remove ${out}_QC4.remove \
      --make-bed --out ${out}_QC5

# check that no sex mismatches and duplicates are left
plink --bfile ${out}_QC5 --extract ${out}_QC2_pruned.prune.in --check-sex 0.3 \
      --out ${out}_QC5_sex
awk '$5 == "PROBLEM"' ${out}_QC5_sex.sexcheck
plink --bfile ${out}_QC5 --extract ${out}_QC2_pruned.prune.in \
      --genome --min 0.125 --out ${out}_QC5_rels_dups
awk '$10 > 0.9' ${out}_QC5_rels_dups.genome



################################################################################
#### Ancestry
################################################################################
#
# Note: Since R ver 4.0, stringsAsFactors is set FALSE by default. This breaks
# the ancestry inference R script in KING because e1071::svm() requires the
# dependent variable to be of type factor and not of type character to correctly
# infer "classification machine" as type. 

# get list of QC filtered and pruned SNPs
plink --bfile ${out}_QC5 --extract ${out}_QC2_pruned.prune.in \
      --make-bed --out ${out}_QC5_pruned
      
# project study sample into PC space of 1000 Genomes reference 
# (--prefix does not work properly if it contains a path, mv output afterwards)
king -b ${king_1kg_ref}.bed,${out}_QC5_pruned.bed \
     --pca --projection --rplot --prefix ${outpref}_QC5_king
mv ${outpref}_QC5_king* ${outdir}
rm .RData

# keep only samples annotated as EUR ancestry to avoid confounding due to
# population stratification
awk 'BEGIN {OFS="\t"} {if ($9 != "EUR") {print $1,$2}}' \
  ${out}_QC5_king_InferredAncestry.txt | tail -n +2 > ${out}_QC5_nonEUR.remove
plink --bfile ${out}_QC5 --remove ${out}_QC5_nonEUR.remove \
      --make-bed --out ${out}_QC6
sed 's/$/\tNon-EUR ancestry/' ${out}_QC5_nonEUR.remove >> ${out}.sample_exclusions.tsv



################################################################################
#### Genetic outliers
################################################################################
# due to necessary temporary removal of relatives, run in two parts: 
#   A) IID1 from each pair of relatives removed
#   B) IID2 from each pair of relatives removed

# generate list of pairs of relatives with pi hat > 0.125
plink --bfile ${out}_QC6 --extract ${out}_QC2_pruned.prune.in \
      --genome --min 0.125 --out ${out}_QC6_rels_dups
awk '{if ($10 <= 0.9) print $1,$2,$3,$4,$10}' ${out}_QC6_rels_dups.genome \
    > ${out}_QC6_rels.flag

cut -f1,2 -d' ' ${out}_QC6_rels.flag | sed 's/ /\t/' > ${out}_QC6_rels.A.temp_remove
cut -f3,4 -d' ' ${out}_QC6_rels.flag | sed 's/ /\t/' > ${out}_QC6_rels.B.temp_remove


# ## Option A: MDS approach like Till Andlauer's v2019 genotype QC
# for part in A B; do
#   plink --bfile ${out}_QC6 --extract ${out}_QC2_pruned.prune.in \
#         --remove ${out}_QC6_rels.${part}.temp_remove \
#         --cluster --mind 0.05 --mds-plot 10 --out ${out}_QC6_outlier.${part}
#   Rscript bin/get_mds_outliers.R 4 \
#       ${out}_QC6_outlier.${part}.mds ${out}_QC6_outlier.${part}
# done
# 
# 
# ## Option B: MDS-based iterative outlier removal
# for part in A B; do
#   touch ${out}_QC6_outlier.${part}.remove 
# 
#   # iterative outlier detection
#   i=1
#   while true; do
#     # combine list of relatives to be removed with detected outliers to be removed
#     cat ${out}_QC6_rels.${part}.temp_remove ${out}_QC6_outlier.${part}.remove |
#       grep -v outlier >> ${out}_QC6_rels_outlier.${part}.remove
#   
#     plink --bfile ${out}_QC6 --extract ${out}_QC2_pruned.prune.in \
#           --remove ${out}_QC6_rels_outlier.${part}.remove \
#           --cluster --mind 0.05 --mds-plot 10 \
#           --out ${out}_QC6_outlier.${part}_round${i}
#     
#     Rscript bin/get_mds_outliers.R 4 \
#         ${out}_QC6_outlier.${part}_round${i}.mds \
#         ${out}_QC6_outlier.${part}_round${i} 
#         
#     if [ -f "${out}_QC6_outlier.${part}_round${i}.remove" ]; then
#       echo "Genetic outlier round $i" >> ${out}_QC6_outlier.${part}.remove
#       cat ${out}_QC6_outlier.${part}_round${i}.remove >> ${out}_QC6_outlier.${part}.remove
#       i=$[$i + 1]
#     else
#       echo "No further outliers found in part $part round $i"
#       break
#     fi
#   done
# done


## Option C: EIGENSOFT approach: iterative removal of genetic outliers, with  
# outlier defined as >6 SD distance from the mean within PC1-10
for part in A B; do
  touch ${out}_QC6_outlier.${part}.remove 

  # iterative outlier detection
  i=1
  while true; do
    # combine list of relatives to be removed with detected outliers to be removed
    cat ${out}_QC6_rels.${part}.temp_remove ${out}_QC6_outlier.${part}.remove |
      grep -v outlier >> ${out}_QC6_rels_outlier.${part}.remove
  
    plink --bfile ${out}_QC6 --extract ${out}_QC2_pruned.prune.in \
          --remove ${out}_QC6_rels_outlier.${part}.remove --pca 10 header tabs \
          --out ${out}_QC6_outlier.${part}_round${i}
    
    Rscript bin/get_pca_outliers.R 6 \
        ${out}_QC6_outlier.${part}_round${i}.eigenvec \
        ${out}_QC6_outlier.${part}_round${i} 
        
    if [ -f "${out}_QC6_outlier.${part}_round${i}.remove" ]; then
      echo "Genetic outlier round $i" >> ${out}_QC6_outlier.${part}.remove
      cat ${out}_QC6_outlier.${part}_round${i}.remove >> ${out}_QC6_outlier.${part}.remove
      i=$[$i + 1]
    else
      echo "No further outliers found in part $part round $i"
      break
    fi
  done
done

cat ${out}_QC6_outlier.{A,B}.remove | grep -v FID | grep -v round | cut -f1,2 |
  sort | uniq > ${out}_QC6_outlier.remove
plink --bfile ${out}_QC6 --remove ${out}_QC6_outlier.remove \
      --make-bed --out ${out}_QC7
sed 's/$/\tGenetic outlier/' ${out}_QC6_outlier.remove >> ${out}.sample_exclusions.tsv



################################################################################
#### Additional variant filtering and ID conversion
################################################################################

# convert Illumina IDs to rs IDs where necessary and exclude duplicate variants
# (-> ${out}_QC7.duplicates.exclude and ${out}_QC7.updatedIDs_dupsFlagged.bim)
Rscript bin/update_Illumina_IDs.R ${out}_QC7.bim $illumina2rsid

# identify ambiguous SNPs
awk '$5=="a" && $6=="t" || 
     $5=="t" && $6=="a" || 
     $5=="c" && $6=="g" || 
     $5=="g" && $6=="c" {print $2}' ${out}_QC7.updatedIDs_dupsFlagged.bim \
     > ${out}_QC7.updatedIDs_dupsFlagged.ambiguous.snps

# combine variant exclusion lists from the last two commands
cat ${out}_QC7.duplicates.exclude \
    ${out}_QC7.updatedIDs_dupsFlagged.ambiguous.snps | sort | uniq \
    > ${out}_QC7.dups_ambiguous.exclude

# remove:
# - duplicate variants (based on identical ID)
# - non-autosomal variants
# - ambiguous SNPs
# - call rate < 98%
# - MAF < 1%
# - HWE p value < 1e-6
plink --bfile ${out}_QC7 --bim ${out}_QC7.updatedIDs_dupsFlagged.bim \
      --chr 1-22 --exclude ${out}_QC7.dups_ambiguous.exclude \
      --geno 0.02 --maf 0.01 --hwe include-nonctrl 1e-6 \
      --make-bed --out ${out}_postQC

# generate final list of pairs of relatives with pi hat > 0.125
plink --bfile ${out}_postQC --extract ${out}_QC2_pruned.prune.in \
      --genome --min 0.125 --out ${out}_postQC_rels

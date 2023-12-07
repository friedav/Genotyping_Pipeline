#!/usr/bin/env bash
#
# Quality control and preprocessing of genotype data from PLINK file set: Part 1
#
# The preimputation genotype QC consists of two parts, since in case of the
# identification of sex mismatches and genetic duplicates in part 1 you might
# want to investigate and resolve some of these issues before continuing with
# part 2, instead of simply excluding all these samples.
#
# Prerequesits: Sex column in FAM file reflects phenotypic sex info if available
#
# Usage:        
# 01_run_genotype_QC_part1.sh \
#       <path to bfile set of input data> \
#       <prefix to be used for processed data> \
#       <output directory>
#
# General thoughts: 
# 1. Is technical data quality okay (variants/individuals)? -> Call rate, Fhet
# 2. Is individual identity as expected? -> Sex mismatches, genetic duplicates
set -e

# mamba env create -f environment.yml -p ./env
conda activate ./env

raw=$1
outpref=$2
outdir=$3

mkdir -p $outdir
out=$outdir/$outpref

# initialize summary of sample exclusions
echo -e "FID\tIID\tReason" > ${out}.sample_exclusions.tsv


################################################################################
#### Pre-filtering: Call rate and MAF
################################################################################

# remove SNPs with low call rate
# (note: order of operations: --missing before --geno/--maf -> *.imiss reflects 
#  individual-level call rates before filtering of variants)
plink --bfile ${raw} --missing --geno 0.02 --maf 0.01 \
      --make-bed --out ${out}_QC1       

# remove individuals with low call rate 
# (note: order of operations: --mind before --missing -> *.lmiss reflects 
#  variant-level call rates after filtering of individuals)
plink --bfile ${out}_QC1 --mind 0.02 --missing --make-bed --out ${out}_QC2
if [ -f "${out}_QC2.irem" ]; then
    sed 's/$/\tCall rate/' ${out}_QC2.irem >> ${out}.sample_exclusions.tsv
fi



################################################################################
#### Pruning
################################################################################

# retrieve list of approximately independent variants (LD pruning) after
# additional variant filtering (including removal of extended MHC region on
# chromosome 6, 25-35 Mb and a typical inversion site on chromosome 8, 7-13 Mb)
# for usage in multiple subsequent steps of this script
plink --bfile ${out}_QC2 --geno 0.02 --hwe include-nonctrl 1e-3 --maf 0.05 \
      --exclude range data/remove_prune.txt \
      --indep-pairphase 200 100 0.2 --out ${out}_QC2_pruned
      


################################################################################
#### Autosomal heterozygosity (Fhet/Inbreeding) - cf. RICOPILI
################################################################################

# Note: if sample is small, take MAF values from reference via --read-freq;
# LD pruning before --het calculation suggested in Plink documentation
plink --bfile ${out}_QC2 --extract ${out}_QC2_pruned.prune.in \
      --het --out ${out}_QC2_het
tail -n +2  ${out}_QC2_het.het |
    awk -v 'OFS=\t' '{if ($6 > 0.2 || $6 < -0.2) print $1,$2}' \
    > ${out}_QC2_het.remove
plink --bfile ${out}_QC2 --remove ${out}_QC2_het.remove \
      --make-bed --out ${out}_QC3
sed 's/$/\tAutosomal heterozygosity/' ${out}_QC2_het.remove \
    >> ${out}.sample_exclusions.tsv



################################################################################
#### Sex 
################################################################################
#
# Note from Plink 1.9 documentation:
# "We suggest running --check-sex once without parameters, eyeballing the 
#  distribution of F estimates (there should be a clear gap between a very tight
#  male clump at the right side of the distribution and the females everywhere 
#  else), and then rerunning with parameters corresponding to the empirical gap."

# identify sex mismatches and samples with X-chromosomal F coefficient <0.8 & >0.3
# (after looking at distribution of F coeffs, increased female thres from 0.2 to
# 0.3, because two samples - depsyumr00641 with F=0.204 and depsyumr01001 with
# F=0.216 - visually still belong to normal distribution of females)
plink --bfile ${out}_QC3 --extract ${out}_QC2_pruned.prune.in --check-sex 0.3 \
      --out ${out}_QC3_sex
awk '$5 == "PROBLEM"' ${out}_QC3_sex.sexcheck

# issues with x chromosomal heterozygosity
tail -n +2  ${out}_QC3_sex.sexcheck | awk -v 'OFS=\t' '{$1=$1};1' |
    awk -v 'OFS=\t' '{if ($6 > 0.3 && $6 < 0.8) print $1,$2}' \
    > ${out}_QC3_xhet.remove
sed 's/$/\tX-chromosomal heterozygosity/' ${out}_QC3_xhet.remove \
    >> ${out}.sample_exclusions.tsv

# remove samples with too high/too low x-chromosomal heterozygosity, but keep 
# sex mismatches for calculation of relatedness, as it might help to resolve 
# genetic duplicates
plink --bfile ${out}_QC3 --remove ${out}_QC3_xhet.remove \
      --make-bed --out ${out}_QC4

# prepare file for removal of sex mismatches after checking genetic duplicates
awk -v 'OFS=\t' '{if ($5 == "PROBLEM") print $1,$2}' ${out}_QC3_sex.sexcheck \
    > ${out}_QC3_sex.remove



################################################################################
#### Genetic relatedness (PI HAT)
################################################################################

# estimate genetic duplicates and genetic relatedness
plink --bfile ${out}_QC4 --extract ${out}_QC2_pruned.prune.in \
      --genome --min 0.125 --out ${out}_QC4_rels_dups
    
# genetic duplicates (PI_HAT > 0.9), due to ...
awk '$10 > 0.9' ${out}_QC4_rels_dups.genome > ${out}_QC4_dups.genome

# ... a) repeated genotyping:
# some duplicates are due to repeated genotyping; these should have the same
# USI and only differ by a suffix (e.g. _re) -> from each set of genetic dups,
# keep only the one sample with best call rate
grep -E '([a-z]{8}[0-9]{5}).+([a-z]{8}[0-9]{5}).+([a-z]{8}[0-9]{5}).+(\2)' \
    ${out}_QC4_dups.genome > ${out}_QC4_dups.same_IID.genome
Rscript bin/compare_callrate.R \
    ${out}_QC2.imiss \
    ${out}_QC4_dups.same_IID.genome \
    ${out}_QC4_dups.notBestCallrate.remove
sed 's/$/\tRepeated genotyping, not best call rate/' \
  ${out}_QC4_dups.notBestCallrate.remove >> ${out}.sample_exclusions.tsv

# ... b) sample swaps:
# these do not have the same USI (find by inverted grep match)
grep -Ev '([a-z]{8}[0-9]{5}).+([a-z]{8}[0-9]{5}).+([a-z]{8}[0-9]{5}).+(\2)' \
     ${out}_QC4_dups.genome | grep -v PI_HAT > ${out}_QC4_dups.different_IID.genome



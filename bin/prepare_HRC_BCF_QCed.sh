#!/bin/bash
#SBATCH --job-name=make_BCF
#SBATCH -p medium
#SBATCH -c 4
#SBATCH -t 1-0
#SBATCH --mem=64G
#SBATCH --array=1-22
#SBATCH --chdir=/home/fdavid/data/HRC1.1
#SBATCH --account=ag_ukihg_forstner
#
# For genotype phasing and imputation, prepare BCF files for autosomes from 
# HRC1.1 reference panel that was provided as VCF files;
# perform additional filtering / quality control of BCFs according to 
# https://dx.doi.org/10.17504/protocols.io.xbgfijw
# (note: not all steps were relevant/applicable; copy-pasted only description 
#  to relevant steps)
#"- Remove the rare variants, here singletons and doubletons by setting AC 
#   threshold with 'bcftools view'.
# - Split multiallelic sites to biallelic records with 'bcftools norm'.
# - Align the variants to reference genome with 'bcftools norm' in order to have
#   the REF and ALT alleles in the shortest possible representation and to 
#   confirm that the REF allele matches the reference genome, additionally 
#   remove duplicate variants (-d none).
# - After alignment, remove multiallelic records with 'bcftools view', since 
#   these are formed during the alignment if the REF does not match with the 
#   reference genome.
# - Finally, remove sites containing missing data with 'bcftools view'."

module load bzip2 

outdir=BCF_QCed
mkdir -p $outdir

fasta=/home/fdavid/data/1kg_Phase3/human_g1k_v37.fasta
chr=$SLURM_ARRAY_TASK_ID
vcf=HRC.r1-1.EGA.GRCh37.chr${chr}.haplotypes*.vcf.gz

bcftools +fill-tags $vcf -Ou  |
    bcftools view -e 'INFO/AC<3 | INFO/AN-INFO/AC<3'  -Ou | \
    bcftools norm -m -any -Ou | \
    bcftools norm -f $fasta -d none -Ou | \
    bcftools view -m 2 -M 2 -Ou | \
    bcftools view -g ^miss -Ob -o ${outdir}/${vcf/.vcf.gz/.bcf}

bcftools index -f ${outdir}/${vcf/.vcf.gz/.bcf}


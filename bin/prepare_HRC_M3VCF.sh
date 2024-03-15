#!/bin/bash
#SBATCH --job-name=make_M3VCF
#SBATCH -p medium
#SBATCH -c 4
#SBATCH -t 1-0
#SBATCH --mem=64G
#SBATCH --array=1-22
#SBATCH --chdir=/home/fdavid/data/HRC1.1
#SBATCH --account=ag_ukihg_forstner
#
# For genotype imputation with Minimac4, prepare M3VCF files for autosomes from 
# HRC1.1 reference panel that was provided as VCF files
#
# Note: 
# "Currently Minimac4 can ONLY handle M3VCF format files. If your reference 
#  panel is in VCF format, please use Minimac3 to convert the VCF file to M3VCF 
#  (along with parameter estimation) and then use that M3VCF for imputation
#  using Minimac4."

chr=$SLURM_ARRAY_TASK_ID
vcf=HRC.r1-1.EGA.GRCh37.chr${chr}.haplotypes*.vcf.gz

minimac3 --refHaps $vcf \
         --processReference \
         --cpus 4 \
         --prefix M3VCF/${vcf/.vcf.gz/}

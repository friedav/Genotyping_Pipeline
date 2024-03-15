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
# For genotype phasing with Eagle2, prepare BCF files for autosomes from HRC1.1 
# reference panel that was provided as VCF files

module load bzip2 

chr=$SLURM_ARRAY_TASK_ID
vcf=HRC.r1-1.EGA.GRCh37.chr${chr}.haplotypes*.vcf.gz

bcftools view -Ob $vcf -o BCF/${vcf/.vcf.gz/.bcf}
bcftools index -f BCF/${vcf/.vcf.gz/.bcf}

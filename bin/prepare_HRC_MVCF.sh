#!/bin/bash
#SBATCH --job-name=make_MVCF
#SBATCH -p medium
#SBATCH -c 4
#SBATCH -t 1-0
#SBATCH --mem=64G
#SBATCH --array=1-22
#SBATCH --chdir=/home/fdavid/data/HRC1.1
#SBATCH --account=ag_ukihg_forstner
#
# For genotype imputation with Minimac4, prepare MVCF files for autosomes from 
# HRC1.1 reference panel that was provided as VCF files

module load CMake/3.23.1-GCCcore-11.3.0  \
            GCC/12.2.0 GCCcore/12.2.0 \
            GLib/2.72.1-GCCcore-11.3.0 \
            libarchive/3.6.1-GCCcore-12.2.0 \
            libffi/3.4.2-GCCcore-11.3.0 \
            zlib/1.2.12-GCCcore-12.2.0

chr=$SLURM_ARRAY_TASK_ID
vcf=HRC.r1-1.EGA.GRCh37.chr${chr}.haplotypes*.vcf.gz

minimac4 --compress-reference $vcf > MVCF/${vcf/.vcf.gz/.msav}

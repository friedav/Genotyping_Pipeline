#!/usr/bin/env bash
# The publicly available part of the Haplotype Reference Consortium (HRC) 
# reference dataset v1.1 was requested via the Wellcome Sanger eDAM portal
# in January 2024 for usage at the Institute of Human Genetics in Bonn
# by Friederike David. Authorised personnal additionally include
# Charlotte Pahnke, Priyatam Dutta, Carina Mathey, and Andreas Forstner.
#
# This file documents the preparation of the data files for use in the
# genotype imputation workflow. Not meant/tested for execution as script.


#### Data retrieval ####

# download dataset from the European Genome-Phenome Archive (EGA) via pyega3
datadir="/home/fdavid/data/HRC1.1"
mkdir $datadir
pyega3 -cf ~/.ega-credentials.json -c 10 fetch EGAD00001002729 --output-dir .

# remove empty tmp_download folders and move all files into one folder
rmdir $datadir/EGAF0000138*/.tmp_download
mv $datadir/EGAF0000138*/* $datadir/
rmdir $datadir/EGAF0000138*

# link to imputation project directory
ln -s $datadir /home/fdavid/projects/Genotyping_Pipeline/data


#### Formatting to bcf ####
# for eagle2; formatting requires bcftools (already installed for imputation.nf)

mkdir ${datadir}/BCF
sbatch /home/fdavid/projects/Genotyping_Pipeline/bin/prepare_HRC_BCF.sh


#### Formatting to MVCF files ####
# for minimac4, according to https://github.com/statgen/Minimac4

cd $datadir
mkdir MVCF
sbatch /home/fdavid/projects/Genotyping_Pipeline/bin/prepare_HRC_MVCF.sh


#### Formatting to M3VCF files [DEPRECATED] ####
# for minimac4

# I overlooked the "IMPORTANT NOTE !!!" on genome.sph.umich.edu/wiki/Minimac4
# stating that this "wiki is outdated and only applicable to version 4.0.x 
# (a.k.a 1.0.x). Documentation for version 4.1.x can be found on the Minimac4 
# Github page", so the following code chunk can be considered deprecated...
#
# # according to https://genome.sph.umich.edu/wiki/Minimac4:
# # "Currently Minimac4 can ONLY handle M3VCF format files. If your reference 
# #  panel is in VCF format, please use Minimac3 to convert the VCF file to M3VCF 
# #  (along with parameter estimation) and then use that M3VCF for imputation
# #  using Minimac4."
# 
# # download link did not work (connection to FTP server not possible)
# #wget ftp://share.sph.umich.edu/minimac3/Minimac3Executable.tar.gz
# 
# # instead, clone from Github and compile from source
# cd ~/bin
# git clone https://github.com/Santy-8128/Minimac3.git
# cd Minimac3
# make
# ln -s /home/fdavid/bin/Minimac3/bin/Minimac3 /home/fdavid/bin/minimac3
# cd -
# 
# # prepare M3VCF files 
# cd $datadir
# mkdir M3VCF
# 
# ## to time consuming, wrote batch script for parallelization
# #for vcf in HRC.r1-1.EGA.GRCh37.chr*.haplotypes*.vcf.gz; do
# #    minimac3 --refHaps $vcf \
# #             --processReference \
# #             --cpus 4 \
# #             --prefix M3VCF/${vcf/.vcf.gz/}
# #done
# 
# sbatch /home/fdavid/projects/Genotyping_Pipeline/bin/prepare_HRC_M3VCF.sh



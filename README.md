# Overview

- Step 0: Manual preparation of raw genotype data as Plink bfile set (*.bed, *.bim, *.fam) with phenotypic sex info in *.fam
- Step 1: Pre-imputation quality control part 1 - [bin/01_run_genotype_QC_part1.sh](bin/01_run_genotype_QC_part1.sh)
- Step 2: Manual preparation of sample exclusions and sample swaps based on sample identity issues detected in part 1 (sex mismatches, duplicates)
- Step 3: Pre-imputation quality control part 2 - [bin/02_run_genotype_QC_part2.sh](bin/02_run_genotype_QC_part2.sh)
- Step 4: Imputation -> [imputation.nf](imputation.nf)

Note: The pre-imputation QC scripts represent templates that should be modified to suit a given dataset / project.


# Imputation

Phasing with Eagle2 and imputation with Minimac4 is implemented in the Nextflow script [imputation.nf](imputation.nf).
As reference data, either the 1000 Genomes Phase 3 (`1kg`) or the Haplotype Reference Consortium (`HRC`) are intended.

Example command to start Nextflow run:

```bash
nextflow run imputation.nf -resume -with-report \
  --outdir "output" \
  --outpref "example" \
  --bfile_dir "output/from/preimputation_QC_part2" \
  --bfile_pref "example_postQC" \
  --ref "1kg"
```

`outdir` and `outpref` specify the output directory and prefix string for the expected output.  
`bfild_dir` and `bfile_pref` specify the directory and prefix of the post-QC pre-imputation Plink bfile set to be used as input.  
`ref` should be "1kg" or "HRC" and specifies which reference panel to use. 


# Requirements / Getting Started

```
git clone https://github.com/friedav/Genotyping_Pipeline.git
cd Genotyping_Pipeline
conda env create -f environment.yml -p ./env
```

Tools that need to be manually installed and present in your `$PATH`:
- Nextflow v22.10.x or earlier (important because script is based on DSL1, which is no longer supported by newer Nextflow versions)
- Eagel v2.4.0 (https://alkesgroup.broadinstitute.org/Eagle/)
- Minimac4 (https://github.com/statgen/Minimac4)
- bcftools (http://samtools.github.io/bcftools/)
- GATK4 (https://github.com/broadinstitute/gatk)
- KING (https://www.kingrelatedness.com/Download.shtml) - only for pre-imputation QC script 2

Tools that need to be manually installed and paths modified in `imputation.nf`:
- `check_script`: Will Rayner's script for "HRC or 1000G Imputation preparation and checking" (https://www.chg.ox.ac.uk/~wrayner/tools/)

Reference data files that need to be manually prepared and paths modified in `imputation.nf`:
- reference haplotypes
  - `refvariants`: variant table, e.g. 1000GP_Phase3_combined.legend based on https://www.chg.ox.ac.uk/~wrayner/tools/1000GP_Phase3_combined.legend.gz
  - `refeagle`: chr-level BCF files to be used by Eagle, e.g. see https://alkesgroup.broadinstitute.org/Eagle/#x1-310005.3.1
  - `refminimac`: chr-level MSAV files to be used by Minimac, e.g. from https://share.sph.umich.edu/minimac4/panels/g1k_p3_msav_files_with_estimates.tar.gz
- `dbsnp`: dbSNP variants for annotation of rsIDs, e.g. from https://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
- `genmap`: genetic coordinates, e.g. those provided with the Eagel installation (Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz)





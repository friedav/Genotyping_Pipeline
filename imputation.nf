#!/usr/bin/env nextflow
/*
 * Genotype imputation
 */
 

/**** Data input and parameters ****/

// output
params.outdir = ""
params.outpref = ""

// use either HRC or 1000 Genomes (1kg) reference data
params.ref = "HRC"

// target post-QC bfile set
params.bfile_dir=""
params.bfile_pref=""
bfile = Channel.value( [params.bfile_pref, 
                       file("${params.bfile_dir}/${params.bfile_pref}.bed"),
                       file("${params.bfile_dir}/${params.bfile_pref}.bim"),
                       file("${params.bfile_dir}/${params.bfile_pref}.fam")] )

if(params.ref =~ /HRC/) {
  // reference haplotypes
  refdir = "$baseDir/data/HRC1.1"
  refvariants = file("$refdir/HRC.r1-1.EGA.GRCh37.sites.tab")
  refeagle = Channel.fromFilePairs("$refdir/BCF/HRC.r1-1.EGA.GRCh37.chr*.haplotypes.{bcf,bcf.csi}")
                    .filter{ it[0] =~ /chr\d+/ }
                    .map{ pref, files -> [(pref =~ /chr(\d+)/)[0][1], files[0], files[1]] }
  refminimac = Channel.fromPath("$refdir/MVCF/HRC.r1-1.EGA.GRCh37.chr*.haplotypes.msav")
                      .filter( ~/.*\d+\.haplotypes.msav/ )
                      .map{ file -> [ (file.name =~ /chr(\d+)/)[0][1], file ] } 

} else if(params.ref =~ /1kg/) {
  // reference haplotypes
  refdir = "$baseDir/data/1kg_Phase3"
  refvariants = file("$refdir/1000GP_Phase3_combined.legend")
  refeagle = Channel.fromFilePairs("$refdir/ALL.chr*.phase3_integrated.20130502.genotypes.{bcf,bcf.csi}")
                    .filter{ it[0] =~ /chr\d+/ }
                    .map{ pref, files -> [(pref =~ /chr(\d+)/)[0][1], files[0], files[1]] }
  refminimac = Channel.fromPath("$refdir/Minimac_G1K_P3_MSAV_FILES_WITH_ESTIMATES/*.1000g.Phase3.v5.With.Parameter.Estimates.msav")
                      .filter( ~/.*\d+\.1000g.Phase3.*/ )
                      .map{ file -> [ (file.name =~ /^\d+/)[0], file ] }    

} else {
  error "Error: Currently, only HRC or 1kg are implemented as --ref parameter"
}

// reference fasta file
reffasta = file("$baseDir/data/1kg_Phase3/human_g1k_v37.fasta")  

// genetic coordinates
genmap = file("/home/fdavid/bin/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz")

// chromosomes to parallelize over (typecast to string necessary to combine with refeagle)
chrs = Channel.of(1..22).map{ it -> "$it" }

// thresholds for post-imputation filtering
params.thres_maf = 0.01
params.thres_r2 = 0.3



/**** Processes ****/

/* prepare target bcf for phasing:
 * split input data by chromosome, recode into vcf, split multiallelic sites, 
 * align to reference, make index
 */
process prep_target_bcf {
  input:
  val chr from chrs  
  tuple val(pref), file(bed), file(bim), file(fam) from bfile
  file reffasta
  
  output:
  tuple val(chr), file("${pref}_chr${chr}.bcf"), file("*.bcf.csi") into target_bcf
  
  script:
  """
  module load bzip2 
 
  # make allele codes uppercase
  awk '{print \$1, \$2, \$3, \$4, toupper(\$5), toupper(\$6)}' $bim \
        > ${pref}.uppercaseAlleles.bim

  plink2 \
    --bfile $pref --chr ${chr}  \
    --bim ${pref}.uppercaseAlleles.bim \
    --snps-only just-acgt \
    --recode vcf id-paste=iid bgz \
    --ref-from-fa $reffasta \
    --out ${pref}_chr${chr}

  # split multiallelic sites into biallelic records
  bcftools view -Ou ${pref}_chr${chr}.vcf.gz | \
    bcftools norm -Ou -m -any | \
    bcftools norm -Ou -d none -f $reffasta --check-ref ws | \
    bcftools view -Ob -o ${pref}_chr${chr}.bcf -m 2 -M 2 && \
    bcftools index -f ${pref}_chr${chr}.bcf
  """
}

  
/* phasing */
process run_phasing {
  tag "chr${chr}"

  input:
  tuple val(chr), file(refbcf), file(refidx), file(targetbcf), file(targetidx) from refeagle.combine(target_bcf, by:0)
  file(genmap)

  output:
  tuple val(chr), file("${params.outpref}_chr${chr}.${params.ref}_phased.vcf.gz") into target_phased

  script:
  """
  eagle --vcfTarget $targetbcf \
        --vcfRef $refbcf \
        --allowRefAltSwap \
        --geneticMapFile $genmap \
        --chrom $chr \
        --outPrefix ${params.outpref}_chr${chr}.${params.ref}_phased \
        --vcfOutFormat z \
        --numThreads ${task.cpus} \
    2>&1 | tee eagle_chr${chr}.log
  """
}


/* imputation */
process run_imputation {
  tag "chr${chr}"
  publishDir "$params.outdir", mode: 'copy'

  input:
  tuple val(chr), file(target), file(refminimac) from target_phased.combine(refminimac, by:0)

  output:
  tuple val(chr), file("${params.outpref}_chr${chr}.${params.ref}_imputed.vcf.gz") into target_imputed
    
  script:
  """
  module load CMake/3.23.1-GCCcore-11.3.0  \
              GCC/12.2.0 GCCcore/12.2.0 \
              GLib/2.72.1-GCCcore-11.3.0 \
              libarchive/3.6.1-GCCcore-12.2.0 \
              libffi/3.4.2-GCCcore-11.3.0 \
              zlib/1.2.12-GCCcore-12.2.0

  # make index file
  bcftools index -t --threads ${task.cpus} $target

  # run imputation
  minimac4 $refminimac $target \
           --format GT,DS,GP \
           --output ${params.outpref}_chr${chr}.${params.ref}_imputed.vcf.gz \
           --output-format vcf.gz \
           --threads ${task.cpus}
  """
}


/* post-imputation filtering by MAF and imputation INFO */
process filter_imputed {
  tag "chr${chr}"
  publishDir "$params.outdir", mode: 'copy'

  input:
  tuple val(chr), file(vcf) from target_imputed

  output:
  file("${params.outpref}_chr${chr}.${params.ref}_imputed.*.vcf.gz") into target_imputed_filtered

  script:
  """
  module load bzip2

  bcftools filter -e "INFO/MAF < $params.thres_maf" -Ou $vcf | \
    bcftools filter -e "INFO/R2 < $params.thres_r2" \
    -Oz -o ${params.outpref}_chr${chr}.${params.ref}_imputed.MAF${params.thres_maf}_R2${params.thres_r2}.vcf.gz
  """
}


/* extract variant info file */
process extract_variants {
  publishDir "$params.outdir", mode: 'copy'
  module 'R'

  input:
  file vcfs from target_imputed_filtered.toList()

  output:
  file("${params.outpref}.${params.ref}_imputed.MAF${params.thres_maf}_R2${params.thres_r2}.variants.tsv")

  script:
  """
  #!/usr/bin/env Rscript
  library(tidyverse); library(data.table)

  # fields in INFO column
  info <- c("TYPE", "AF", "MAF", "AVG_CS", "R2", "ER2")

  variants <- list.files(pattern = ".vcf.gz") %>%
    lapply(function(file) {
      df <- fread(cmd = paste('zgrep -v "##"', file, '| cut -f1-8')) %>%
        mutate(INFO = str_replace(INFO, "TYPED;IMPUTED", "TYPED"))
  }) %>% bind_rows() %>%
    separate(INFO, into = info, sep = ";") %>%
    mutate(across(any_of(info[-1]), ~str_remove(.x, "^.+=") %>% as.numeric())) %>%
    arrange(`#CHROM`, POS)

  fwrite(variants, "${params.outpref}.${params.ref}_imputed.MAF${params.thres_maf}_R2${params.thres_r2}.variants.tsv", sep = "\\t")
  """
}


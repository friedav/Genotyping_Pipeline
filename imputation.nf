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
  // flag for William Rayner's perl script whether HRC or 1000G is used
  check_script_refflag = "--hrc"

  // reference haplotypes
  refdir = "$baseDir/data/HRC1.1"
  refvariants = file("$refdir/HRC.r1-1.EGA.GRCh37.sites.tab")
  refeagle = Channel.fromFilePairs("$refdir/BCF_QCed/HRC.r1-1.EGA.GRCh37.chr*.haplotypes.{bcf,bcf.csi}")
                    .filter{ it[0] =~ /chr\d+/ }
                    .map{ pref, files -> [(pref =~ /chr(\d+)/)[0][1], files[0], files[1]] }
  refminimac = Channel.fromPath("$refdir/MVCF/HRC.r1-1.EGA.GRCh37.chr*.haplotypes.msav")
                      .filter( ~/.*\d+\.haplotypes.msav/ )
                      .map{ file -> [ (file.name =~ /chr(\d+)/)[0][1], file ] } 

} else if(params.ref =~ /1kg/) {
  // flag for William Rayner's perl script whether HRC or 1000G is used
  check_script_refflag = "--1000g"

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

// genetic coordinates
genmap = file("/home/fdavid/bin/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz")

// William Rayner's perl script to align target to reference data
check_script = "$baseDir/bin/HRC-1000G-check/HRC-1000G-check-bim.pl"

// thresholds for post-imputation filtering
params.thres_maf = 0.01
params.thres_r2 = 0.3



/**** Processes ****/

/* run pre-imputation checks on target and reference data using Will Rayner's
 * script from https://www.chg.ox.ac.uk/~wrayner/tools/, according to
 * https://imputationserver.readthedocs.io/en/latest/prepare-your-data/
 */
process check_target {
  input:
  tuple val(pref), file(bed), file(bim), file(fam) from bfile
  file refvariants

  output:
  file("*check-updated*.vcf") into target_checked

  script:
  """
  # make allele codes uppercase
  awk '{print \$1, \$2, \$3, \$4, toupper(\$5), toupper(\$6)}' $bim \
        > ${pref}.uppercaseAlleles.bim

  plink \
    --bfile $pref \
    --bim ${pref}.uppercaseAlleles.bim \
    --freq \
    --make-bed \
    --out ${pref}_check

  perl $check_script -b ${pref}_check.bim \
                     -f ${pref}_check.frq \
                     -r $refvariants "$check_script_refflag"

  sh ./Run-plink.sh
  """
}

target_checked.flatten()
              .map{ file -> [ (file.name =~ /chr(\d+).vcf/)[0][1], file ] }
              .set{ target_checked_byChr }


/* prepare target bcf for phasing and make bcf index
 */ 
process prep_target_bcf {
  tag "chr${chr}"

  input:
  tuple val(chr), file(vcf) from target_checked_byChr
    
  output:
  tuple val(chr), file("*.bcf"), file("*.bcf.csi") into target_bcf
  
  script:
  """
  module load bzip2 

  # convert to bcf (keep only biallelic sites) and make index
  bcftools view -Ob $vcf -o ${params.outpref}_chr${chr}.bcf -m 2 -M 2
  bcftools index -f ${params.outpref}_chr${chr}.bcf
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
           --all-typed-sites \
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
  tuple val(chr), 
    file("${params.outpref}_chr${chr}.${params.ref}_imputed.*.vcf.gz") into target_imputed_filtered

  script:
  """
  module load bzip2

  bcftools filter -e "INFO/MAF < $params.thres_maf" -Ou $vcf | \
    bcftools filter -e "INFO/R2 < $params.thres_r2" \
    -Oz -o ${params.outpref}_chr${chr}.${params.ref}_imputed.MAF${params.thres_maf}_R2${params.thres_r2}.vcf.gz
  """
}


/* extract variant info files */
process extract_variants {
  input:
  tuple val(chr), file(vcf) from target_imputed_filtered

  output:
  file("variants*.fixedQUAL.tsv") into variants_by_chr

  script:
  """
  module load Java/17.0.6
  module load bzip2

  bcftools view --drop-genotypes -Oz -o variants_chr${chr}.vcf.gz $vcf
  bcftools index --tbi -f variants_chr${chr}.vcf.gz

  gatk VariantsToTable -V variants_chr${chr}.vcf.gz  -O variants_chr${chr}.tsv

  # for some unknown reason, GATK's VariantsToTable replaces the "." in QUAL 
  # with "-10.0"; use awk to revert back to expected "."
  # (see https://github.com/broadinstitute/gatk/issues/8748)
  awk -F "\\t" 'BEGIN {OFS=FS} {
    if ( \$6 == "-10.0" ) \$6 = "."
    print \$0
    }' variants_chr${chr}.tsv > variants_chr${chr}.fixedQUAL.tsv
  """
}


/* merge chr-level variant info files */
process join_variant_tables {
  publishDir "$params.outdir", mode: 'copy'

  input:
  file(variants) from variants_by_chr.toSortedList( { a, b -> 
                            def chr_a = (a =~ /chr(\d+)/)[0][1] as Integer
                            def chr_b = (b =~ /chr(\d+)/)[0][1] as Integer
                            return chr_a <=> chr_b } )

  output:
  file("${params.outpref}.${params.ref}_imputed.MAF${params.thres_maf}_R2${params.thres_r2}.variants.tsv")

  script:
  """
  awk 'NR==1 || FNR > 1' $variants > \
    "${params.outpref}.${params.ref}_imputed.MAF${params.thres_maf}_R2${params.thres_r2}.variants.tsv"
  """
}

#!/usr/bin/env nextflow

params.help = false
params.reads = null
params.ref_dir = "./"
params.seqtype = "dna"
params.singleEnd = false
params.cpu = 4
params..vcf_dir = "./vcf"


/* Prints help when asked for and exits */
def helpMessage() {
    log.info"""
    =========================================
    neoflow => WXS anylsis
    =========================================
    Usage:
    nextflow run neoflow_.vcf.nf
    Arguments:
      --reads                     Reads data in fastq.gz or fastq format. For example, "*_{1,2}.fastq.gz"
      --ref_dir                   HLA reference folder
      --seqtype                   Read type, dna or rna. Default is dna.
      --singleEnd                 Single end or not, default is false (pair end reads)
      --cpu                       The number of CPUs, default is 4.
      --vcf_dir                   Folder of variant file , default is "./"
      --help                      Print help message

    """.stripIndent()
}
// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

/* Click  vcf_dir folder has been created  */
vcf_dir = file(params.vcf_dir)
if(!vcf_dir.isDirectory()){
	vcf_dir_result = vcf_dir.mkdirs()
	println vcf_dir_result ? "Create folder: $vcf_dir!" : "Cannot create directory: $myDir!"
}



reference_dir = file(params.ref_dir)

/* Click pair-end or single-end sequence and input sequence data  */
Channel.fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
	.ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs" + "to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { input_data }

/* Sequence alignment by bwa */
process reads_mapping{

        echo true

        tag "reads_mapping"
        cpus "$params.cpu"

        input:
        file reference_dir
        set val(pattern), file(reads) from input_data

        output:
        set val(pattern), "${pattern}.sam" into mapped_reads

        script:

        """
        bwa mem -p -t ${params.cpu} -M ${reference_dir}/GRCh37.fa \
        -R "@RG\\tID:sample_1\\tLB:sample_1\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:sample_1" \
         ${reads[0]} ${reads[1]}  > ${pattern}.sam

        """
}



/* Duplicates Marking by Picard Tools
 * (only run one time)*/

process run_make_duplicates{
        cpus "$params.cpu"
        input:
        file reference_dir
        set val(pattern), "${pattern}.sam"  from mapped_reads

        output:
        set val(pattern), "${pattern}_dup_sort.bam" into reads_mapping_duplicates_sort_1,reads_mapping_duplicates_sort_2
        set val(pattern), "${pattern}_depth_out.txt" into reads_mapping_duplicates_depth_out

        script:

        """
        samtools view -@ ${params.cpu}  -S -b ${pattern}.sam > ${pattern}.bam
        samtools sort -n  -@ ${params.cpu} -o ${pattern}.sort.bam ${pattern}.bam
        samtools fixmate -@ ${params.cpu} -m  ${pattern}.sort.bam ${pattern}.fixmate.bam
        samtools sort   -@ ${params.cpu} -o ${pattern}.positionsort.bam  ${pattern}.fixmate.bam
        samtools markdup  ${pattern}.positionsort.bam ${pattern}_dup_sort.bam
        samtools index ${pattern}_dup_sort.bam
        samtools depth -a ${pattern}_dup_sort.bam >  ${pattern}_depth_out.txt
        """
}

/*  Call Variants by gatk4 (before BQSR) */

process call_variant{
        cpus "$params.cpu"
        input:
        file reference_dir
        set val(pattern), "${pattern}_dup_sort.bam" from reads_mapping_duplicates_sort_1

        output:
        set val(pattern), "${pattern}_raw_variants.vcf" into raw_variants

        script:
        """
        gatk    HaplotypeCaller  \
        -R ${reference_dir}/GRCh37.fa \
        -I ${pattern}_dup_sort.bam \
        -O ${pattern}_raw_variants.vcf
        """
}


/* Extract SNPs and Indels vcf by gatk4 (before BQSR) */

process select_variant{

        input:
        file reference_dir
        set val(pattern), "${pattern}_raw_variants.vcf" from raw_variants
        output:
        set val(pattern), "${pattern}_raw_snps.vcf" into raw_variants_snps
        set val(pattern), "${pattern}_raw_indels.vcf" into raw_variants_indels

        script:
        """
        gatk  SelectVariants \
        -R ${reference_dir}/GRCh37.fa \
        -V ${pattern}$_raw_variants.vcf \
        --select-type-to-include  SNP \
        -O ${pattern}_raw_snps.vcf
        gatk SelectVariants \
        -R ${reference_dir}/GRCh37.fa \
        -V ${pattern}_raw_variants.vcf \
        --select-type-to-include  INDEL \
        -O ${pattern}_raw_indels.vcf
        """
}

 /* Filter_variants by gatk4 (before BQSR) */

 process filter_variants{

         input:
         file reference_dir
         set val(pattern), "${pattern}_raw_snps.vcf" from raw_variants_snps
         set val(pattern), "${pattern}_raw_indels.vcf" from  raw_variants_indels

         output:
         set val(pattern), "${pattern}_filtered_snps.vcf" into filtered_variants_snps
         set val(pattern), "${pattern}_filtered_indels.vcf" into  filtered_variants_indels

         script:
         """
         #snp
         gatk VariantFiltration \
         -R ${reference_dir}/GRCh37.fa \
         -V ${pattern}_raw_snps.vcf \
         --filter-name "QD_filter" \
         --filter-expression "QD < 2.0" \
         --filter-name "FS_filter" \
         --filter-expression "FS > 60.0" \
         --filter-name  "MQ_filter" \
         --filter-expression "MQ < 40.0" \
         --filter-name "SOR_filter" \
         --filter-expression "SOR > 10.0" \
         -O ${pattern}_filtered_snps.vcf

         #Indels
         gatk  VariantFiltration \
         -R ${reference_dir}/GRCh37.fa \
         -V ${pattern}_raw_indels.vcf \
         --filter-name "QD_filter" \
         --filter-expression "QD < 2.0" \
         --filter-name "FS_filter" \
         --filter-expression "FS > 60.0" \
         --filter-name "MQ_filter" \
         --filter-expression "MQ < 40.0" \
         --filter-name "SOR_filter" \
         --filter-expression "SOR > 10.0" \
         -O ${pattern}_filtered_indels.vcf
         """
 }

/* Exclude Filtered Variants(get BQSR vcf) */

process get_bqsr_.vcf{

        input:
        file reference_dir
        set val(pattern), "${pattern}_filtered_snps.vcf" from filtered_variants_snps
        set val(pattern), "${pattern}_filtered_indels.vcf"  from filtered_variants_indels

        output:
        set val(pattern), "${pattern}_bqsr_snps.vcf" into bqsr_variants_snps
        set val(pattern), "${pattern}_bqsr_indels.vcf"  into bqsr_variants_indels

        script:
        """
        gatk SelectVariants \
        -R ${reference_dir}/GRCh37.fa \
        -V ${pattern}_filtered_snps.vcf \
        --select-type-to-include SNP \
        -O ${pattern}_bqsr_snps.vcf

        gatk SelectVariants \
        -R ${reference_dir}/GRCh37.fa \
        -V ${pattern}_filtered_indels.vcf \
        --select-type-to-include INDEL \
        -O ${pattern}_bqsr_indels.vcf
        """
}

/*Base Quality Score Recalibration (BQSR)*/
process BQSR{

        input:
        file reference_dir
        set val(pattern), "${pattern}_bqsr_snps.vcf" from bqsr_variants_snps
        set val(pattern), "${pattern}_bqsr_indels.vcf"  from  bqsr_variants_indels
        set val(pattern), "${pattern}_dup_sort.bam" from reads_mapping_duplicates_sort_2
        output:
        set val(pattern), "${pattern}_recal_reads.bam" into bqsr_read_map
        script:
        """
        gatk BaseRecalibrator \
        -R ${reference_dir}/GRCh37.fa \
        -I ${pattern}_dup_sort.bam \
        --known-sites ${pattern}_bqsr_snps.vcf \
        --known-sites ${pattern}_bqsr_indels.vcf \
        -O ${pattern}_recal_data.table

        gatk ApplyBQSR \
        -R ${reference_dir}/GRCh37.fa \
        -I ${pattern}_dup_sort.bam \
        -bqsr ${pattern}_recal_data.table \
        -O ${pattern}_recal_reads.bam \
        """
}

/*  Call Variants by gatk4 (after BQSR) */

process call_variant_bqsr{
        cpus "$params.cpu"
        input:
        file reference_dir
        set val(pattern), "${pattern}_recal_reads.bam" from bqsr_read_map

        output:
        set val(pattern), "${pattern}_raw_variants_recal.vcf" into variants_recal
        script:
        """
        gatk   HaplotypeCaller  \
        -R ${reference_dir}/GRCh37.fa \
        -I ${pattern}_recal_reads.bam \
        -O ${pattern}_raw_variants_recal.vcf
        """
}

/* Extract SNPs and Indels vcf by gatk4 (after BQSR) */

process select_variant_bqsr{

        input:
        file reference_dir
        set val(pattern), "${pattern}._raw_variants_recal.vcf" from variants_recal

        output:
        set val(pattern), "${pattern}_raw_snps_recal.vcf" into variants_snps_recal
        set val(pattern), "${pattern}_raw_indels_recal.vcf" into variants_indels_recal

        script:
        """
        gatk  SelectVariants \
        -R ${reference_dir}/GRCh37.fa \
        -V ${pattern}_raw_variants_recal.vcf \
        --select-type-to-include SNP \
        -O ${pattern}_raw_snps_recal.vcf
        gatk  SelectVariants \
        -R ${reference_dir}/GRCh37.fa \
        -V ${pattern}_raw_variants_recal.g.vcf \
        --select-type-to-include SNP \
        -O ${pattern}_raw_snps_recal.vcf
        """
}

/* Filter_variants by gatk4 (after BQSR) */

process filter_variants_bqsr{
        publishDir "${vcf_dir}/${pattern}/", mode: "copy", overwrite: true
        input:
        file reference_dir
        set val(pattern), "${pattern}_raw_snps_recal.vcf" from variants_snps_recal
        set val(pattern), "${pattern}_raw_indels_recal.vcf" from variants_indels_recal

        output:
        set val(pattern), "${pattern}_snps.vcf" into variants_snps

        script:
        """
        #snp
        gatk VariantFiltration \
        -R ${reference_dir}/GRCh37.fa \
        -V ${pattern}_raw_snps_recal.vcf \
        --filter-name "QD_filter" \
        --filter-expression "QD < 2.0" \
        --filter-name "FS_filter" \
        --filter-expression "FS > 60.0" \
        --filter-name "MQ_filter" \
        --filter-expression "MQ < 40.0" \
        --filter-name "SOR_filter" \
        --filter-expression "SOR > 10.0" \
        -O ${pattern}_snp.vcf
        """
}

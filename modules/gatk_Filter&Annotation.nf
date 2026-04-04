#!/usr/bin/env nextflow

process gatk_extract_filtered_variants {
    
    publishDir "${params.variant_call}/filtered_variants", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (cohort_metadata), path (filtered_vcf), path (filtered_vcf_index)
    output:
    tuple val (cohort_metadata), path ("${filtered_vcf.simpleName}_extracted.vcf.gz"), emit: "extracted_filtered_variants"
    tuple val (cohort_metadata), path ("${filtered_vcf.simpleName}_extracted.vcf.gz.tbi"), emit: "extracted_filtered_variants_index"

    script:

    """
    gatk SelectVariants \
    -R ${params.ref} \
    -V ${filtered_vcf} \
    --exclude-filtered \
    -O ${filtered_vcf.simpleName}_extracted.vcf.gz
    """

}



process DOWNLOAD_SNPEFF_DB {
    storeDir "${params.snpeff_db_dir}/snpeff_cache" // This ensures the download is saved permanently
    conda "bioconda::snpeff=5.4.0c"

    output:
    path "${params.snpeff_db}", emit: snpeff_db_path

    script:
    """
    snpEff download -v ${params.snpeff_db} -dataDir \$(pwd)
    """
}

process SNPEFF_ANNOTATE {
    tag "${filtered_extracted_vcf.simpleName}"
    publishDir "${params.outdir}/annotation", mode: 'copy'
    conda "bioconda::snpeff=5.4.0c bioconda::htslib=1.19"

    input:
    tuple val(cohort_metadata), path(filtered_extracted_vcf)
    path(snpeff_db) // This comes from the download process

    output:
    tuple val(cohort_metadata), path("${filtered_extracted_vcf.simpleName}_ann.vcf.gz"), emit: ann_vcf
    tuple val(cohort_metadata), path("${filtered_extracted_vcf.simpleName}_ann.vcf.gz.tbi"), emit: ann_vcf_index

    script:
    """
    # -dataDir . tells snpEff to look in the current working directory 
    # (where Nextflow linked the db_dir)
    snpEff \
    -Xmx8g \
    -dataDir \$(pwd) \
    ${params.snpeff_db} \
    ${filtered_extracted_vcf} | bgzip > ${filtered_extracted_vcf.simpleName}_ann.vcf.gz
    tabix -p vcf ${filtered_extracted_vcf.simpleName}_ann.vcf.gz
    """
}

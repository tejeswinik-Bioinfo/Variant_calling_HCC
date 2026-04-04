#!/usr/bin/env nextflow

process gatk_HaplotypeCaller{

    publishDir "${params.outdir}/variants", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (metadata), path (bqsr_bam)
    tuple val (metadata), path (bqsr_bam_index)
    output:
    tuple val (metadata), path ("*.g.vcf.gz"), emit: "gvcf"
    tuple val (metadata), path ("*.g.vcf.gz.tbi"), emit: "gvcf_index"

    script:
    sample_id = metadata.sampleName
    """
    gatk HaplotypeCaller \
    -R ${params.ref} \
    -I ${bqsr_bam} \
    -ERC GVCF \
    -O ${sample_id}.g.vcf.gz 
    """
}

process gatk_combine_gvcfs {

    publishDir "${params.outdir}/variants", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    path (gvcf_file)
    path (gvcf_index)

    output:
    path ("*.g.vcf.gz"), emit: "vcf"
    path ("*.g.vcf.gz.tbi"), emit: "vcf_index"

    script:
    
    def variants_list = gvcf_file.collect{ "--variant ${it}" }.join(' ')
    """
    gatk CombineGVCFs \
    -R ${params.ref} \
    ${variants_list} \
    -O ${params.cohort}_combined.g.vcf.gz

    """
}

process gatk_GenotypeGVCFs {

    publishDir "${params.outdir}/variants", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (cohort_metadata), path (gvcf_file)
    tuple val (cohort_metadata), path (gvcf_index)

    output:
    tuple val (cohort_metadata), path ("*.vcf.gz"), emit: "vcf"
    tuple val (cohort_metadata), path ("*.vcf.gz.tbi"), emit: "vcf_index"

    script:
    
    """
    gatk GenotypeGVCFs \
    -R ${params.ref} \
    -V ${gvcf_file} \
    -O ${params.cohort}_joint.vcf.gz

    """
}

process gatk_select_variants_SNPs {

    publishDir "${params.variant_call}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (cohort_metadata), path (jointgvcf_file)
    tuple val (cohort_metadata), path (jointgvcf_file_index)
    output:
    tuple val (cohort_metadata), path ("*_snp.vcf.gz"), emit: "snp_vcf"
    tuple val (cohort_metadata), path ("*_snp.vcf.gz.tbi"), emit: "snp_vcf_index"

    script:
    cohort_id = cohort_metadata.cohortName
    """
    gatk SelectVariants \
    -R ${params.ref} \
    -V ${jointgvcf_file} \
    --select-type-to-include SNP \
    -O ${cohort_id}_snp.vcf.gz
    """

}

process gatk_select_variants_INDELs {

    publishDir "${params.variant_call}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (cohort_metadata), path (jointgvcf_file)
    tuple val (cohort_metadata), path (jointgvcf_index)
    output:
    tuple val (cohort_metadata), path ("*_indels.vcf.gz"), emit: "indel_vcf"
    tuple val (cohort_metadata), path ("*_indels.vcf.gz.tbi"), emit: "indel_vcf_index"

    script:
    
    """
    gatk SelectVariants \
    -R ${params.ref} \
    -V ${jointgvcf_file} \
    --select-type-to-include INDEL \
    -O ${params.cohort}_indels.vcf.gz
    """

}

process gatk_FilterVariants_SNPs {

    publishDir "${params.variant_call}/filtered_variants", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (cohort_metadata), path(jointgvcf_file)
    tuple val (cohort_metadata), path(jointgvcf_index)
    output:
    tuple val (cohort_metadata), path ("*_filtered_snps.vcf.gz"), emit: "filtered_snp"
    tuple val (cohort_metadata), path ("*_filtered_snps.vcf.gz.tbi"), emit: "filtered_snp_index"

    script:
    
    """
    gatk VariantFiltration \
    -R ${params.ref} \
    -V ${jointgvcf_file} \
    --filter-expression "QD < 2.0 || QUAL < 30.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    -filter-name 'myfilter' \
    -O ${params.cohort}_filtered_snps.vcf.gz 
    """
}

process gatk_FilterVariants_INDELs {
    publishDir "${params.variant_call}/filtered_variants", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (cohort_metadata), path (indel_vcf_file)
    tuple val (cohort_metadata), path (indel_vcf_index)
    output:
    tuple val (cohort_metadata), path("*_filtered_indels.vcf.gz"), emit: "filtered_indels"
    tuple val (cohort_metadata), path("*_filtered_indels.vcf.gz.tbi"), emit: "filtered_indels_index"

    script:

    """
    gatk VariantFiltration \
    -R ${params.ref} \
    -V ${indel_vcf_file} \
    --filter-expression "QD < 2.0 || QUAL < 30.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
    -filter-name 'myfilter' \
    -O ${params.cohort}_filtered_indels.vcf.gz
    """
}




process gatk_extract_filtered_indels {
    publishDir "${params.variant_call}/filtered_variants", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (cohort_metadata), path (filtered_indels_vcf) 
    tuple val (cohort_metadata), path (filtered_indels_vcf_index)
    output:
    tuple val (cohort_metadata), path ("*_extracted_filtered_indels.vcf.gz"), emit: "extracted_filtered_indels"
    tuple val (cohort_metadata), path ("*_extracted_filtered_indels.vcf.gz.tbi"), emit: "extracted_filtered_indels_index"

    script:

    """
    gatk SelectVariants \
    -R ${params.ref} \
    -V ${filtered_indels_vcf} \
    --exclude-filtered \
    -select 
    """
}
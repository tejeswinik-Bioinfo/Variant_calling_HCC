#!/usr/bin/env nextflow

process gatk_mark_duplicates {

    publishDir "${params.outdir}/aligned_reads", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (metadata), path (aligned_bam)
    output:
    tuple val (metadata), path ("*_sorted_dedup*.bam"), emit: "dedup_bam"
    tuple val (metadata), path ("*_sorted_dedup*.bai"), emit: "dedup_index"

    script:
    sample_id = metadata.sampleName
    """
    gatk MarkDuplicatesSpark \
    -I ${aligned_bam} \
    -O ${sample_id}_sorted_dedup_reads.bam
    """
}

process gatk_base_recalibrator {

    publishDir "${params.outdir}/aligned_reads", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val(metadata), path (dedup_bam)
    output:
    tuple val (metadata), path ("*_recal_data.table")

    script:
    sample_id = metadata.sampleName
    """
    gatk BaseRecalibrator \
    -I ${dedup_bam} \
    -R ${params.ref} \
    --known-sites ${params.known_sites} \
    -O ${sample_id}_recal_data.table
    """
}

process gatk_applybqsr{

    publishDir "${params.outdir}/aligned_reads", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (metadata), path (dedup_bam), path (recal_table)
    output:
    tuple val (metadata), path ("*dedup_bqsr*.bam"), emit: "dedup_bqsr_bam"
    tuple val (metadata), path ("*dedup_bqsr*.bai"), emit: "dedup_bqsr_index"

    script:
    sample_id = metadata.sampleName
    """
    gatk ApplyBQSR \
    -I ${dedup_bam} \
    -R ${params.ref} \
    --bqsr-recal-file ${recal_table} \
    -O ${sample_id}_dedup_bqsr.bam
    """

}








#!/usr/bin/env nextflow

process gatk_mark_duplicates {

    publishDir "${params.outdir}/aligned_reads", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (metadata), path (aligned_bam)
    output:
    tuple val (metadata), path ("*_sorted_dedup*")

    script:
    sample_id = metadata.sampleName
    """
    gatk MarkDuplicates \
    -I ${aligned_bam} \
    -O ${sample_id}_sorted_dedup_reads.bam
    """
}





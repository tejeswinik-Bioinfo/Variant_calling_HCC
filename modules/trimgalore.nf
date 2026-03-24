#!/usr/bin/env nextflow

process trimgalore_trim_reads{
    publishDir "${params.trimmed_fastq}/trimmed_fastq", mode: "copy"
    conda "bioconda::trim-galore=0.6.11"

    input:
    tuple val(metadata), path (r1), path (r2)
    output:
    tuple val(metadata), path ("*trimmed_R1.fastq.gz"), path ("*trimmed_R2.fastq.gz"), emit: "trimmed_fastq"
    path ("*.html"), emit: "trimmed_reports"

    script:
    sample_id = metadata.sampleName
    """
    trim_galore \
    --paired \
    --quality 30 \
    --length 25 \
    --cores ${params.threads} \
    "${r1}" "${r2}"
    """
}
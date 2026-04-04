#!/usr/bin/env nextflow

process bwa_index{

    storeDir "${params.ref_parent}"
    conda "bioconda::bwa=0.7.19"

    input:
    path (ref)
    output:
    path ("${ref}.*")

    script:
    ref_parent = file(params.ref).getParent()

    """
    bwa index ${ref}
    """
}


process bwa_align_reads{

    publishDir "${params.outdir}/aligned_reads", mode: "copy"
    conda "bioconda::bwa=0.7.19 bioconda::samtools=1.23.1"

    input:
    tuple val(metadata), path (r1), path (r2)
    output:
    tuple val(metadata), path ("*_aligned_reads*")

    script:
    sample_id = metadata.sampleName
    """
    bwa mem \
    -t ${params.threads} \
    -R "@RG\\tID:${sample_id}\\tPL:ILUMINA\\tSM:${sample_id}" \
    ${params.ref} \
    ${r1} \
    ${r2} | samtools view -bS - > ${sample_id}_aligned_reads.bam
    """    
}
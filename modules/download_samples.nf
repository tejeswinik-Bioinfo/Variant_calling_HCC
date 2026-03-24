#!/usr/bin/env nextflow

process prefetch_download_samples{

    publishDir "/mnt/data/biouser/raw_data/", mode: "copy"
    conda "bioconda::sra-tools=3.2.1"

    input:
    tuple val(metadata), val(sample_id)

    output:
    tuple val(metadata), path("${sample_id}.sra"), emit: sra_files

    script:
    """
    prefetch --option-file ${sample_id} -X u
    """
}

process fasterq_dump_DownloadAndSplit{
    publishDir "/mnt/data/biouser/raw_data/HCC_fastq/", mode: "copy"
    conda "bioconda::sra-tools=3.2.1"

    input:
    tuple val(metadata), val(sra_file)

    output:
    tuple val(metadata), path("${sra_file}_R1.fastq.gz"), path("${sra_file}_R2.fastq.gz")

    script:
    sample_id = metadata.sampleName

    """
    fasterq-dump \
    --split-files ${sra_file} \
    --progress

    mv ${sample_id}_1.fastq ${sample_id}_R1.fastq
    mv ${sample_id}_2.fastq ${sample_id}_R2.fastq
    
    gzip *.fastq
    """

}





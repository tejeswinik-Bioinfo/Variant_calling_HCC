#!/usr/bin/env nextflow

params.input = ""
params.mode = ""
params.sampleName = ""
params.cohort = ""
params.r1 = ""
params.r2 = ""


if (params.mode == 'single') {

    if (!params.sampleName) {
    error "Parameter --sampleName is required for single-sample mode."
}
        if (!params.r1 || !params.r2) {
            error "Single-sample mode requires --r1 <path> and --r2 <path>."
    }
            def meta = [sampleName: params.sampleName, pairedEnd: true]
            sample_ch = Channel.of([meta, file(params.r1), file(params.r2)])



} else if (params.mode == 'multi') {
    if (!params.cohort) {
        error "Parameter --cohort is required for multi-sample mode."
    }
        if (!params.input) {
            error "Multi-sample mode requires --input <samplesheet.csv>."
        }
        sample_ch = Channel.fromPath(params.input)
            .splitCsv(header: true, sep: ",")
            .map { row ->
                def meta = [sampleName: row.sample_name, pairedEnd: row.paired]
                def r1 = row.file_path_R1
                def r2 = row.file_path_R2
                return [meta, r1, r2]
            }
          
} else {
    error "Unknown params.mode value '${params.mode}'. Use 'single' or 'multi'."
}


include {fasterq_dump_DownloadAndSplit} from "./modules/download_samples.nf"
include {trimgalore_trim_reads} from "./modules/trimgalore"
include {bwa_index; bwa_align_reads} from "./modules/bwa"
include {gatk_mark_duplicates; gatk_base_recalibrator; gatk_applybqsr} from "./modules/gatk_preprocessing"
include {gatk_HaplotypeCaller; gatk_combine_gvcfs; gatk_GenotypeGVCFs; gatk_select_variants_SNPs; gatk_select_variants_INDELs; gatk_FilterVariants_SNPs; gatk_FilterVariants_INDELs} from "./modules/gatk_variant_calling"
include {gatk_extract_filtered_variants;DOWNLOAD_SNPEFF_DB;SNPEFF_ANNOTATE} from "./modules/gatk_Filter&Annotation"


workflow {
    if (params.mode == 'single') {
        single_sample()
    } else {
        multi_samples()
    }
}


workflow common_preprocessing {
    take:
    sample_ch

    main:
    trimmed_reads = trimgalore_trim_reads(sample_ch)
    bwa_index(params.ref)
    bwa_align_reads(trimgalore_trim_reads.out.trimmed_fastq)
    gatk_mark_duplicates(bwa_align_reads.out)
    gatk_base_recalibrator(gatk_mark_duplicates.out.dedup_bam)

    applybqsr_input = gatk_mark_duplicates.out.dedup_bam.join(gatk_base_recalibrator.out)
    gatk_applybqsr(applybqsr_input)
    gatk_HaplotypeCaller(gatk_applybqsr.out.dedup_bqsr_bam, gatk_applybqsr.out.dedup_bqsr_index)

    emit:
    gvcf = gatk_HaplotypeCaller.out.gvcf
    gvcf_index = gatk_HaplotypeCaller.out.gvcf_index
}


workflow single_sample {
    main:
    gvcf_output = common_preprocessing(sample_ch)

    gvcf_ch = gvcf_output.gvcf
    gvcf_index_ch = gvcf_output.gvcf_index

    gatk_GenotypeGVCFs(gvcf_ch, gvcf_index_ch)
    gatk_select_variants_SNPs(gatk_GenotypeGVCFs.out.vcf, gatk_GenotypeGVCFs.out.vcf_index)
    gatk_select_variants_INDELs(gatk_GenotypeGVCFs.out.vcf, gatk_GenotypeGVCFs.out.vcf_index)
    gatk_FilterVariants_SNPs(gatk_select_variants_SNPs.out.snp_vcf, gatk_select_variants_SNPs.out.snp_vcf_index)
    gatk_FilterVariants_INDELs(gatk_select_variants_INDELs.out.indel_vcf, gatk_select_variants_INDELs.out.indel_vcf_index)

    snp_ch = gatk_FilterVariants_SNPs.out.filtered_snp.join(gatk_FilterVariants_SNPs.out.filtered_snp_index)
    indel_ch = gatk_FilterVariants_INDELs.out.filtered_indels.join(gatk_FilterVariants_INDELs.out.filtered_indels_index)

    filter_ch = snp_ch.concat(indel_ch)
    gatk_extract_filtered_variants(filter_ch)

    DOWNLOAD_SNPEFF_DB()
    SNPEFF_ANNOTATE(gatk_extract_filtered_variants.out.extracted_filtered_variants, DOWNLOAD_SNPEFF_DB.out.snpeff_db_path)

    emit:
    annotated_variants = SNPEFF_ANNOTATE.out.ann_vcf
}


workflow multi_samples {
    main:
    gvcf_output = common_preprocessing(sample_ch)

    vcf_list = gvcf_output.gvcf.map { it[1] }.collect()
    tbi_list = gvcf_output.gvcf_index.map { it[1] }.collect()

    gatk_combine_gvcfs(vcf_list, tbi_list)

    cohort_metadata = [cohortName: "${params.cohort}"]
    vcf_meta = gatk_combine_gvcfs.out.vcf.map { vcf -> [cohort_metadata, vcf] }
    tbi_meta = gatk_combine_gvcfs.out.vcf_index.map { tbi -> [cohort_metadata, tbi] }

    gatk_GenotypeGVCFs(vcf_meta, tbi_meta)
    gatk_select_variants_SNPs(gatk_GenotypeGVCFs.out.vcf, gatk_GenotypeGVCFs.out.vcf_index)
    gatk_select_variants_INDELs(gatk_GenotypeGVCFs.out.vcf, gatk_GenotypeGVCFs.out.vcf_index)
    gatk_FilterVariants_SNPs(gatk_select_variants_SNPs.out.snp_vcf, gatk_select_variants_SNPs.out.snp_vcf_index)
    gatk_FilterVariants_INDELs(gatk_select_variants_INDELs.out.indel_vcf, gatk_select_variants_INDELs.out.indel_vcf_index)

    snp_ch = gatk_FilterVariants_SNPs.out.filtered_snp.join(gatk_FilterVariants_SNPs.out.filtered_snp_index)
    indel_ch = gatk_FilterVariants_INDELs.out.filtered_indels.join(gatk_FilterVariants_INDELs.out.filtered_indels_index)

    filter_ch = snp_ch.concat(indel_ch)
    gatk_extract_filtered_variants(filter_ch)

    DOWNLOAD_SNPEFF_DB()
    SNPEFF_ANNOTATE(gatk_extract_filtered_variants.out.extracted_filtered_variants, DOWNLOAD_SNPEFF_DB.out.snpeff_db_path)

    emit:
    annotated_variants = SNPEFF_ANNOTATE.out.ann_vcf
}

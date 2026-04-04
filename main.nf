#!/usr/bin/env nextflow

params.input = ""
include {fasterq_dump_DownloadAndSplit} from "./modules/download_samples.nf"
include{trimgalore_trim_reads} from "./modules/trimgalore"
include{bwa_index; bwa_align_reads} from "./modules/bwa"
include{gatk_mark_duplicates; gatk_base_recalibrator; gatk_applybqsr} from "./modules/gatk_preprocessing"
include{gatk_HaplotypeCaller; gatk_combine_gvcfs; gatk_GenotypeGVCFs; gatk_select_variants_SNPs; gatk_select_variants_INDELs; gatk_FilterVariants_SNPs; gatk_FilterVariants_INDELs} from "./modules/gatk_variant_calling"
include{gatk_extract_filtered_variants;DOWNLOAD_SNPEFF_DB;SNPEFF_ANNOTATE} from "./modules/gatk_Filter&Annotation"


workflow {  

    sample_ch = Channel.fromPath(params.input)
    .splitCsv(header: true, sep: ",")
    .map{row ->
        def meta = [sampleName: row.sample_name, pairedEnd: row.paired]
        def r1 = row.file_path_R1
        def r2 = row.file_path_R2
        return [meta,r1,r2]

        }


    //fasterq_dump_DownloadAndSplit(sample_ch).view() 
    trimmed_reads = trimgalore_trim_reads(sample_ch)
    


    //bwa_index_files = [file("${params.ref}.ann"), file("${params.ref}.amb"), file("${params.ref}.bwt"), file("${params.ref}.pac"), file("${params.ref}.sa")]
    //index_exists = bwa_index_files.every{it.exists()}

    //if (!index_exists) {
    //    bwa_index(params.ref)
    //}
    //else {
        //println "BWA index files exixt. Skipping BWA indexing step."
    //}
    

    bwa_index(params.ref)

    bwa_align_reads(trimgalore_trim_reads.out.trimmed_fastq)

    gatk_mark_duplicates(bwa_align_reads.out)
  
    gatk_base_recalibrator(gatk_mark_duplicates.out.dedup_bam)

    applybqsr_input = gatk_mark_duplicates.out.dedup_bam.join(gatk_base_recalibrator.out)
    gatk_applybqsr(applybqsr_input)

    gatk_HaplotypeCaller(gatk_applybqsr.out.dedup_bqsr_bam, gatk_applybqsr.out.dedup_bqsr_index)

   //Collect gvcf files and their indexes for CombineGVCFs
    vcf_list = gatk_HaplotypeCaller.out.gvcf.map { it[1] }.collect()
    tbi_list = gatk_HaplotypeCaller.out.gvcf_index.map { it[1] }.collect()

    gatk_combine_gvcfs(vcf_list,tbi_list)

// Adding cohortName to the metadata for downstream processes
    cohort_metadata = [cohortName: "${params.cohort}"]

    vcf_meta = gatk_combine_gvcfs.out.vcf.map {vcf -> [cohort_metadata,vcf] }
    tbi_meta = gatk_combine_gvcfs.out.vcf_index.map {tbi -> [cohort_metadata,tbi] }
    
    gatk_GenotypeGVCFs(vcf_meta,tbi_meta)

    gatk_select_variants_SNPs(gatk_GenotypeGVCFs.out.vcf, gatk_GenotypeGVCFs.out.vcf_index)
    
    gatk_select_variants_INDELs(gatk_GenotypeGVCFs.out.vcf, gatk_GenotypeGVCFs.out.vcf_index)
    
    gatk_FilterVariants_SNPs(gatk_select_variants_SNPs.out.snp_vcf, gatk_select_variants_SNPs.out.snp_vcf_index)
    gatk_FilterVariants_INDELs(gatk_select_variants_INDELs.out.indel_vcf, gatk_select_variants_INDELs.out.indel_vcf_index)
    
    snp_ch = gatk_FilterVariants_SNPs.out.filtered_snp.join(gatk_FilterVariants_SNPs.out.filtered_snp_index)
    indel_ch = gatk_FilterVariants_INDELs.out.filtered_indels.join(gatk_FilterVariants_INDELs.out.filtered_indels_index)
   
    filter_ch = snp_ch.concat(indel_ch)
    gatk_extract_filtered_variants(filter_ch)
    gatk_extract_filtered_variants.out.extracted_filtered_variants.view()

     // 1. Download DB once
    DOWNLOAD_SNPEFF_DB()
    
    // 2. Run annotation using the output of the download
    SNPEFF_ANNOTATE(gatk_extract_filtered_variants.out.extracted_filtered_variants, DOWNLOAD_SNPEFF_DB.out.snpeff_db_path)

    emit:
    annotated_variants = SNPEFF_ANNOTATE.out.ann_vcf

}

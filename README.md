Introduction
Description
This workflow focuses on Germline short variant discovery. It uses paired-end FASTQ files as input. It consists of pipelines from GATK best practices data pre-processing and Germline short variant discovery. 
Execution Modes:
Single sample
Multi-sample (cohort)

Workflow Architecture
2.1 Pipeline Steps
Quality Control (FastQC)
Trimming (if applicable)
Alignment (BWA-MEM)
Sorting & MarkDuplicates
Base Quality Score Recalibration (BQSR)
Variant Calling (HaplotypeCaller)
Joint Genotyping (for multi-sample)
Variant Filtering
Annotation
Directory Structure
├── main.nf
├── modules
│   ├── bwa.nf
│   ├── download_samples.nf
│   ├── gatk_Filter&Annotation.nf
│   ├── gatk_preprocessing.nf
│   ├── gatk_variant_calling.nf
│   └── trimgalore.nf
├── nextflow.config
├── README.md
└── scripts
    └── ann_plots.py

Requirements
Java 17 or later up to 26
Tools used
Trimgalore
BWA
GATK
SnpEff
MultiQC
Input Requirements
Single sample
Paired-end FASTQ Format:
sample1_R1.fastq.gz
sample1_R2.fastq.gz
Sample Sheet (for multi-sample)
sample_id
paired
file_path_read_1
file_path_read_2
S1
True
/file/path/S1_R1
/file/path/S1_R2
S2
True
/file/path/S2_R1
/file/path/S2_R2


Configuration
7.1 Parameters (params.config)
params {
    outdir = "results/"
    ref = "data/reference/human_genome/genome.fa"
    ref_parent =  "data/reference/”
    known_sites = "data/reference/dbsnp.vcf"
    variant_call = "results/variants"
    cohort = "HCC”
    snpeff_db = "GRCh38.115"
    snpeff_db_dir = "/data/reference/SnpEff"
    run_type = "single" // or "multi”
    threads = 8
}
7.2 Profiles
profiles {
    conda {
        conda.enabled = true
    }
}

7. Running the Pipeline
7.1 Basic Command
nextflow run Variant_calling_HCC --input raw_data/samplesheet_1.csv -profile conda -resume
7.2 Germline – Single Sample
nextflow run main.nf \
  --mode germline \
  --run_type single \
  --input "data/sample/*"
7.3 Germline – Multi Sample
nextflow run main.nf \
  --mode germline \
  --run_type multi \
  --samplesheet samples.csv

8. Output Structure
├── aligned_reads
├── annotation
│   ├── plots
│   └── stats
│       └── plots
├── plots
├── trimmed_reads
└── variants
    └── filtered_variants



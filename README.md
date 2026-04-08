# Germline Variant Calling Workflow

## 1. Introduction

### Description
This workflow focuses on **germline short variant discovery** using **paired-end FASTQ files** as input.  
It implements pipelines based on **GATK Best Practices**, including:

- Data preprocessing  
- Germline variant discovery  

---

## 2. Execution Modes

- **Single Sample**
- **Multi-Sample (Cohort)**

---

## 3. Workflow Architecture

### 3.1 Pipeline Steps

1. Quality Control (**FastQC**)  
2. Trimming (if applicable)  
3. Alignment (**BWA-MEM**)  
4. Sorting & MarkDuplicates  
5. Base Quality Score Recalibration (**BQSR**)  
6. Variant Calling (**HaplotypeCaller**)  
7. Joint Genotyping (for multi-sample)  
8. Variant Filtering  
9. Annotation  

---

## 4. Directory Structure
```text
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

---

## 5. Requirements

- **Java 17 to 26**
---

### 6. Tools Used

- TrimGalore  
- BWA  
- GATK  
- SnpEff  
- MultiQC  

---

## 7. Input Requirements

### 7.1 Single Sample

Paired-end FASTQ files:
sample1_R1.fastq.gz
sample1_R2.fastq.gz


---

### 7.2 Multi-Sample (Sample Sheet)

| sample_id | paired | file_path_read_1 | file_path_read_2 |
|----------|--------|------------------|------------------|
| S1       | True   | /file/path/S1_R1 | /file/path/S1_R2 |
| S2       | True   | /file/path/S2_R1 | /file/path/S2_R2 |

---

## 8. Configuration

### 8.1 Parameters (`params.config`)

```groovy
params {
    outdir = "results/"
    ref = "data/reference/human_genome/genome.fa"
    ref_parent = "data/reference/"
    known_sites = "data/reference/dbsnp.vcf"
    variant_call = "results/variants"
    cohort = "HCC"
    snpeff_db = "GRCh38.115"
    snpeff_db_dir = "/data/reference/SnpEff"
    run_type = "single" // or "multi"
    threads = 8
}

### 8.2 Profiles
profiles {
    conda {
        conda.enabled = true
    }
}

---

## 9. Running the Pipeline
### 9.1 Basic Command
nextflow run Variant_calling_HCC \
  --input raw_data/samplesheet_1.csv \
  -profile conda \
  -resume

### 9.2 Germline – Single Sample
nextflow run main.nf \
  --mode germline \
  --run_type single \
  --input "data/sample/*"
### 9.3 Germline – Multi-Sample
nextflow run main.nf \
  --mode germline \
  --run_type multi \
  --samplesheet samples.csv

## 10. Output Structure

```text
  results/
├── aligned_reads/      # BQSR-calibrated BAMs
├── annotation/         # SnpEff VCFs & Reports
│   ├── plots/          # Variant distribution plots
│   └── stats/          # Annotation statistics
├── plots/              # MultiQC/FastQC reports
├── trimmed_reads/      # Cleaned FASTQ files
└── variants/           # Final VCF files
    └── filtered_variants/

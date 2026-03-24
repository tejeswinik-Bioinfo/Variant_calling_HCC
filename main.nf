#!/usr/bin/env nextflow

params.input = ""
include {fasterq_dump_DownloadAndSplit} from "./modules/download_samples.nf"
include{trimgalore_trim_reads} from "./modules/trimgalore"





workflow {

    sample_ch = Channel.fromPath(params.input)
    .splitCsv(header: true, sep: ",")
    .map{row ->
        def meta = [sampleName: row.sample_name, pairedEnd: row.paired]
        def r1 = row.file_path_R1
        def r2 = row.file_path_R2
        return [meta,r1,r2]

        }
    sample_ch.view()

    //fasterq_dump_DownloadAndSplit(sample_ch).view()
    trimmed_reads = trimgalore_trim_reads(sample_ch)
    trimmed_reads.trimmed_fastq.view()


}
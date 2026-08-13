
include { coveragePlot } from '../modules/generatePlots.nf'
include { getMappedReads } from '../modules/getMappedReads.nf'
include { getUnmappedReads } from '../modules/getUnmappedReads.nf'

workflow GENPLOTS {
    take:
        bams_ch // meta, bam_file, bai_file, is_paired_end
    main:
    // Create sub-channels for each process type
    coverage_ch = bams_ch.map { meta, bam, bai, _is_pe -> tuple(meta, bam, bai) }
    reads_ch = bams_ch.map { meta, bam, _bai, is_pe -> tuple(meta, bam, is_pe) }

    //QC
    //Rendering the depth coverage plot
    coveragePlot(coverage_ch)
    // Check if there are mapped reads
    coveragePlot_out_ch = coveragePlot.out.result
    coveragePlot_out_ch
    | view { result -> log.warn("${result.text}") }

    if ((params.writeMappedReads == true)){
        // write mapped reads
        getMappedReads(reads_ch)

        // write unmappped reads
        getUnmappedReads(reads_ch)
    }
}

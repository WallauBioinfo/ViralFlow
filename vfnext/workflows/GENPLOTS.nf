
include { coveragePlot } from '../modules/generatePlots.nf'
include { getMappedReads } from '../modules/getMappedReads.nf'
include { getUnmappedReads } from '../modules/getUnmappedReads.nf'

workflow GENPLOTS {
    take: 
        bams_ch // sample_id, bam_files, bai_files, is_paired_end
    main:
    //QC
    //Rendering the depth coverage plot
    coveragePlot(bams_ch)
    // Check if there are mapped reads
    coveragePlot_out_ch = coveragePlot.out.result
    coveragePlot_out_ch
    | view( log.warn("${it.text}"))

    if ((params.writeMappedReads == true)){
        // write mapped reads
        getMappedReads(bams_ch)
  
        // write unmappped reads
        getUnmappedReads(bams_ch)
    }
}
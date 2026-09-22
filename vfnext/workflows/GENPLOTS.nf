
include { coveragePlot } from '../modules/generatePlots.nf'
include { getMappedReads } from '../modules/getMappedReads.nf'
include { getUnmappedReads } from '../modules/getUnmappedReads.nf'

workflow GENPLOTS {
    take:
        bamsCh // meta, bam_file, bai_file, is_paired_end
    main:
    // Create sub-channels for each process type
    coverageCh = bamsCh.map { meta, bam, bai, _is_pe -> tuple(meta, bam, bai) }
    readsCh = bamsCh.map { meta, bam, _bai, is_pe -> tuple(meta, bam, is_pe) }

    //QC
    //Rendering the depth coverage plot
    coveragePlot(coverageCh)
    // Check if there are mapped reads
    coveragePlotOutCh = coveragePlot.out.result
    coveragePlotOutCh
    | view { result -> log.warn("${result.text}") }

    if ((params.writeMappedReads == true)){
        // write mapped reads
        getMappedReads(readsCh)

        // write unmappped reads
        getUnmappedReads(readsCh)
    }
}


include { coveragePlot } from '../modules/generatePlots.nf'
include { getMappedReads } from '../modules/getMappedReads.nf'
include { getUnmappedReads } from '../modules/getUnmappedReads.nf'
include { writeMappedReadsEnabled } from '../modules/param_helpers.nf'

workflow GENPLOTS {
    take:
        bamsCh // meta, bam_file, bai_file, is_paired_end
        // The mode's consensus depth threshold, for the coverage plot's stats:
        // params.depth for ILLUMINA, params.np_min_depth for NANOPORE. Taken
        // from the caller rather than read here, because the two modes name it
        // differently and plotting ILLUMINA's `depth` over a NANOPORE run
        // reported recovery at a threshold its consensus never used.
        coverageThreshold
    main:
    // Create sub-channels for each process type
    coverageCh = bamsCh.map { meta, bam, bai, _is_pe -> tuple(meta, bam, bai) }
    readsCh = bamsCh.map { meta, bam, _bai, is_pe -> tuple(meta, bam, is_pe) }

    //QC
    //Rendering the depth coverage plot
    coveragePlot(coverageCh, coverageThreshold)
    // Check if there are mapped reads
    coveragePlotOutCh = coveragePlot.out.result
    coveragePlotOutCh
    | view { result -> log.warn("${result.text}") }

    if (writeMappedReadsEnabled(params.writeMappedReads)) {
        // write mapped reads
        getMappedReads(readsCh)

        // write unmappped reads
        getUnmappedReads(readsCh)
    }
}

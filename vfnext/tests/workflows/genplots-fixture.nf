nextflow.enable.dsl = 2

include { getMappedReads } from '../../modules/getMappedReads.nf'
include { getUnmappedReads } from '../../modules/getUnmappedReads.nf'
include { coveragePlot } from '../../modules/generatePlots.nf'

// Deliberately named something other than "<id>.sorted.bam". Both processes
// used to rebuild that name from meta.id instead of reading the staged path,
// so any BAM the NANOPORE workflow rebound - "<id>.trim.sorted.bam" from
// bamUtil, "<id>.primer_clip.bam" from ampliconclip - made them fail with
// "No such file or directory".
process prepareRenamedBam {
    label "NP_basecontainer"
    tag "${meta.id}"

    input:
        tuple val(meta), path(sam)

    output:
        tuple val(meta),
              path("${meta.id}.trim.sorted.bam"),
              path("${meta.id}.trim.sorted.bam.bai")

    script:
    """
    set -euo pipefail

    samtools view -bS ${sam} \
        | samtools sort -o ${meta.id}.trim.sorted.bam

    samtools index ${meta.id}.trim.sorted.bam
    """
}

workflow GENPLOTS_FIXTURE {
    take:
        samCh

    main:
        prepareRenamedBam(samCh)

        prepareRenamedBam.out
            .map { meta, bam, _bai -> tuple(meta, bam, false) }
            .set { readsCh }

        getMappedReads(readsCh)
        getUnmappedReads(readsCh)
        coveragePlot(prepareRenamedBam.out, params.fixture_coverage_threshold)

    emit:
        mapped = getMappedReads.out
        unmapped = getUnmappedReads.out
        plots = coveragePlot.out.plots
        result = coveragePlot.out.result
}

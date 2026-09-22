nextflow.enable.dsl = 2

include { getMappedReads } from '../../modules/getMappedReads.nf'
include { getUnmappedReads } from '../../modules/getUnmappedReads.nf'

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
        tuple val(meta), path("${meta.id}.trim.sorted.bam")

    script:
    """
    set -euo pipefail

    samtools view -bS ${sam} \
        | samtools sort -o ${meta.id}.trim.sorted.bam
    """
}

workflow GENPLOTS_FIXTURE {
    take:
        samCh

    main:
        prepareRenamedBam(samCh)

        prepareRenamedBam.out
            .map { meta, bam -> tuple(meta, bam, false) }
            .set { readsCh }

        getMappedReads(readsCh)
        getUnmappedReads(readsCh)

    emit:
        mapped = getMappedReads.out
        unmapped = getUnmappedReads.out
}

nextflow.enable.dsl = 2

include { runAmpliconClip } from '../../modules/runAmpliconClip.nf'

process prepareFixtureBam {
    label "NP_basecontainer"
    tag "${meta.id}"

    input:
        tuple val(meta), path(sam)

    output:
        tuple val(meta),
              path("${meta.id}.sorted.bam"),
              path("${meta.id}.sorted.bam.bai")

    script:
    """
    set -euo pipefail

    samtools view -bS ${sam} \
        | samtools sort -o ${meta.id}.sorted.bam

    samtools index ${meta.id}.sorted.bam
    """
}

// nf-test cannot read BAM, so render the clipped alignment as text for the
// assertions.
process dumpSam {
    label "NP_basecontainer"
    tag "${meta.id}"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
        path("${meta.id}.clipped.sam")

    script:
    """
    set -euo pipefail

    samtools view ${bam} > ${meta.id}.clipped.sam
    """
}

workflow AMPLICONCLIP_FIXTURE {
    take:
        samCh
        primer_bed

    main:
        prepareFixtureBam(samCh)
        runAmpliconClip(prepareFixtureBam.out, primer_bed)
        dumpSam(runAmpliconClip.out.bams)

    emit:
        clipped = dumpSam.out
        stats = runAmpliconClip.out.stats
}

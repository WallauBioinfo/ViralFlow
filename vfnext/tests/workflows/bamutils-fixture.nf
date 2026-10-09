nextflow.enable.dsl = 2

include { runBamUtils } from '../../modules/runBamUtils.nf'

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

// nf-test cannot read BAM, so render the trimmed alignment as text.
process dumpSam {
    label "NP_basecontainer"
    tag "${meta.id}"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
        path("${meta.id}.trimmed.sam")

    script:
    """
    set -euo pipefail

    samtools view ${bam} > ${meta.id}.trimmed.sam
    """
}

workflow BAMUTILS_FIXTURE {
    take:
        samCh
        trim_len

    main:
        prepareFixtureBam(samCh)
        runBamUtils(prepareFixtureBam.out, trim_len)
        dumpSam(runBamUtils.out.bams)

    emit:
        trimmed = dumpSam.out
}

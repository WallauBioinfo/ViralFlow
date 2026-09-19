nextflow.enable.dsl = 2

include { run_bam_utils } from '../../modules/runBamUtils.nf'

process prepare_fixture_bam {
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
process dump_sam {
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
        sam_ch
        trim_len

    main:
        prepare_fixture_bam(sam_ch)
        run_bam_utils(prepare_fixture_bam.out, trim_len)
        dump_sam(run_bam_utils.out.bams)

    emit:
        trimmed = dump_sam.out
}

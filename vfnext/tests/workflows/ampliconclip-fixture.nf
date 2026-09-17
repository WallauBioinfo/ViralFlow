nextflow.enable.dsl = 2

include { run_amplicon_clip } from '../../modules/runAmpliconClip.nf'

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

// nf-test cannot read BAM, so render the clipped alignment as text for the
// assertions.
process dump_sam {
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
        sam_ch
        primer_bed

    main:
        prepare_fixture_bam(sam_ch)
        run_amplicon_clip(prepare_fixture_bam.out, primer_bed)
        dump_sam(run_amplicon_clip.out.bams)

    emit:
        clipped = dump_sam.out
        stats = run_amplicon_clip.out.stats
}

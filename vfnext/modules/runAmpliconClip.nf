
process run_amplicon_clip {
    label "NP_basecontainer"
    publishDir { "${params.outDir}/${meta.id}_results/" }, mode: 'copy', overwrite: true
    tag "${meta.id}"

    input:
        tuple val(meta), path(bam), path(bai)
        path(primer_bed)

    output:
        tuple val(meta),
              path("${meta.id}.primer_clip.bam"),
              path("${meta.id}.primer_clip.bam.bai"), emit: bams
        path("${meta.id}.ampliconclip.txt"), emit: stats

    script:
    """
    set -euo pipefail

    # ampliconclip emits reads in an arbitrary order, so sort and index before
    # anything downstream reads depth or calls variants from this BAM.
    samtools ampliconclip \
        --strand \
        --hard-clip \
        -b ${primer_bed} \
        -f ${meta.id}.ampliconclip.txt \
        ${bam} \
        | samtools sort -o ${meta.id}.primer_clip.bam

    samtools index ${meta.id}.primer_clip.bam
    """
}

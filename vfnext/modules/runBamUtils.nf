
process run_bam_utils {
    label "NP_basecontainer"
    publishDir { "${params.outDir}/${meta.id}_results/" }, mode: 'copy', overwrite: true
    tag "${meta.id}"

    input:
        tuple val(meta), path(bam), path(bai)
        val(trim_len)

    output:
        tuple val(meta),
              path("${meta.id}.trim.sorted.bam"),
              path("${meta.id}.trim.sorted.bam.bai"), emit: bams

    script:
    """
    set -euo pipefail

    # --clip soft-clips the trimmed bases. bamUtil's default instead masks them
    # (bases to N, qualities to !) while the alignment keeps spanning the same
    # reference positions, which had two consequences: the masked bases still
    # counted toward the depth consensus masking is computed from while
    # supporting no allele, and Clair3's pileup aborted on them outright:
    #
    #   munmap_chunk(): invalid pointer     (SIGABRT, deterministic)
    #
    # Soft clipping moves the alignment start past the trimmed bases, so they
    # no longer occupy reference positions at all. A trim that would consume a
    # whole read leaves it unclipped and marked unmapped rather than producing
    # a degenerate alignment.
    #
    # Reads are single-end here, so the mate fields bamUtil leaves untouched
    # when clipping need no samtools fixmate pass.
    bam trimBam ${bam} ${meta.id}.trim.bam -L ${trim_len} -R ${trim_len} --clip

    # Required rather than defensive: soft clipping shifts start positions, so
    # trimBam's output is no longer coordinate-sorted, and downstream depth and
    # variant calling both need a sorted, indexed BAM.
    samtools sort ${meta.id}.trim.bam -o ${meta.id}.trim.sorted.bam
    samtools index ${meta.id}.trim.sorted.bam
    """
}

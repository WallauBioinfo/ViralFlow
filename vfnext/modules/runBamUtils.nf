
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

    # bamUtil's default trimBam masks the trimmed bases (sets them to N with
    # quality !) rather than soft-clipping them, so the alignment length is
    # unchanged and the masked bases simply stop supporting any allele.
    bam trimBam ${bam} ${meta.id}.trim.bam -L ${trim_len} -R ${trim_len}

    # Re-sort defensively: trimBam preserves input order, but downstream depth
    # and variant calling both require a coordinate-sorted, indexed BAM.
    samtools sort ${meta.id}.trim.bam -o ${meta.id}.trim.sorted.bam
    samtools index ${meta.id}.trim.sorted.bam
    """
}

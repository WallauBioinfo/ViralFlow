process getUnmappedReads {
    tag "${meta.id}"
    publishDir { "${params.outDir}/${meta.id}_results/" }, mode: "copy"
    label "singlethread"
    input:
        tuple val(meta), path(bam), val(is_paired_end)

    output:
        tuple val(meta), path("*.unmapped.*.fq.gz")

    script:
    """
    # pipefail: without it, a samtools sort that dies partway through its
    # output leaves samtools fastq exiting 0 on the part it got - see
    # getMappedReads.nf, where that was reproduced.
    set -euo pipefail

    samtools view -b -f 4 ${bam} > unmapped.bam
    if [[ ${is_paired_end}  == true ]]; then
      samtools sort -n unmapped.bam | \
      samtools fastq -f 4 -1 ${meta.id}.unmapped.R1.fq.gz -2 ${meta.id}.unmapped.R2.fq.gz
    else
      samtools sort -n unmapped.bam | \
      samtools fastq -f 4 -0 ${meta.id}.unmapped.SE.fq.gz
    fi
    """
}

/*
// --- DOCUMENTATION ----
The mirror of getMappedReads, keeping the reads that did not align.

The single-end branch uses "-0" for the same reason it does there: "-s" is for
singletons, which single-end reads are not, so it would publish an empty
archive and lose the reads to stdout. See getMappedReads.nf for the detail and
for the paired-branch limitation this process shares.
*/

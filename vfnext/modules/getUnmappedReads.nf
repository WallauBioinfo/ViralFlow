process getUnmappedReads {
    publishDir "${params.outDir}/${meta.id}_results/", mode: "copy"
    label "singlethread"
    input:
        tuple val(meta), path(bam_files), val(is_paired_end)

    output:
        tuple val(meta), path("*.unmapped.*.fq.gz")

    script:
    """
    samtools view -b -f 4 ${meta.id}.sorted.bam > unmapped.bam
    if [[ ${is_paired_end}  == true ]]; then
      samtools sort -n unmapped.bam | \
      samtools fastq -f 4 -1 ${meta.id}.unmapped.R1.fq.gz -2 ${meta.id}.unmapped.R2.fq.gz
    else
      samtools sort -n unmapped.bam | \
      samtools fastq -f 4 -s ${meta.id}.unmapped.SE.fq.gz
    fi
    """
}

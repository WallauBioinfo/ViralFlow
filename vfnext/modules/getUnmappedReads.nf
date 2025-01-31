process getUnmappedReads {
    publishDir "${params.outDir}/${sample_id}_results/", mode: "copy"
    label "singlethread"
    input:
        tuple val(sample_id), path(bam_files), val(is_paired_end)

    output:
        tuple val(sample_id), path("*.unmapped.*.fq.gz")

    script:
    """
    samtools view -b -f 4 ${sample_id}.sorted.bam > unmapped.bam
    if [[ ${is_paired_end}  == true ]]; then
      samtools sort -n unmapped.bam | \
      samtools fastq -f 4 -1 ${sample_id}.unmapped.R1.fq.gz -2 ${sample_id}.unmapped.R2.fq.gz
    else
      samtools sort -n unmapped.bam | \
      samtools fastq -f 4 -s ${sample_id}.unmapped.SE.fq.gz
    fi
    """
}

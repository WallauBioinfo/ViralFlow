process getMappedReads{
  tag "${meta.id}"
  publishDir "${params.outDir}/${meta.id}_results/", mode: "copy"
  label "singlethread"
  input:
    tuple val(meta), path(bam_files), val(is_paired_end)
  
  output:
    tuple val(meta), path("*.mapped.*.fq.gz")
  script:
    """
    if [[ ${is_paired_end}  == true ]]; then
      samtools sort -n ${meta.id}.sorted.bam | \
      samtools fastq -F 4 -1 ${meta.id}.mapped.R1.fq.gz -2 ${meta.id}.mapped.R2.fq.gz
    else
      samtools sort -n ${meta.id}.sorted.bam | \
      samtools fastq -F 4 > ${meta.id}.mapped.SE.fq
      gzip ${meta.id}.mapped.SE.fq
    fi
    """
}

/*
// --- DOCUMENTATION ----
This process was designed to get fastqs containing only the mapped reads

1 - "samtools sort -n" organize reads by name (the "sorted.bam" is organized by quality)
2 - "samtools fastq" filter the mapped reads and write it as fastq 
*/

process getMappedReads{
  publishDir "${params.outDir}/${sample_id}_results/", mode: "copy"
  label "singlethread"
  input:
    tuple val(sample_id), path(bam_files), val(is_paired_end)
  
  output:
    tuple val(sample_id), path("*.mapped.*.fq.gz")
  script:
    """
    if [[ ${is_paired_end}  == true ]]; then
      samtools sort -n ${sample_id}.sorted.bam | \
      samtools fastq -F 4 -1 ${sample_id}.mapped.R1.fq.gz -2 ${sample_id}.mapped.R2.fq.gz
    else
      samtools sort -n ${sample_id}.sorted.bam | \
      samtools fastq -F 4 -s ${sample_id}.mapped.SE.fq.gz 
    fi
    """
}

/*
// --- DOCUMENTATION ----
This process was designed to get fastqs containing only the mapped reads

1 - "samtools sort -n" organize reads by name (the "sorted.bam" is organized by quality)
2 - "samtools fastq" filter the mapped reads and write it as fastq 
*/

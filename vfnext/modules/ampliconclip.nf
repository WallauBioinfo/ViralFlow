
process runAmpliconClipping {
  publishDir "${params.outDir}/${sample_id}_results/", mode: "copy"
  
  input:
    tuple val(sample_id), path(bam), path(bai), val(is_paired_end)

  output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), val(is_paired_end), emit: regular_output
    path("${sample_id}.trimmed_reads.txt"), emit: trimmed_reads
    
  script:
    trim_bam = "${sample_id}.trimmed"
    bed = "${params.primersBED}"

    """
    samtools ampliconclip --both-ends --hard-clip \
      --filter-len ${params.minLen} \
      -b ${bed} ${bam} \
      -f ${sample_id}.trimmed_reads.txt > ${trim_bam}

    samtools sort -o ${trim_bam}.sorted.bam ${trim_bam}
    samtools index ${trim_bam}.sorted.bam

    mv ${bam} ${sample_id}.raw.sorted.bam
    mv ${bai} ${sample_id}.raw.sorted.bam.bai
    mv ${trim_bam}.sorted.bam ${sample_id}.sorted.bam
    mv ${trim_bam}.sorted.bam.bai ${sample_id}.sorted.bam.bai
    """
}

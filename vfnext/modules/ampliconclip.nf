
process runAmpliconClipping {
  publishDir "${params.outDir}/${meta.id}_results/", mode: "copy"
  
  input:
    tuple val(meta), path(bam), path(bai), val(is_paired_end)
    path(primers_bed)
    val(min_len)

  output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), val(is_paired_end), emit: regular_output
    path("${meta.id}.trimmed_reads.txt"), emit: trimmed_reads
    
  script:
    sample_id = meta.id
    trim_bam = "${sample_id}.trimmed"
    bed = "${primers_bed}"

    """
    samtools ampliconclip --both-ends --hard-clip \
      --filter-len ${min_len} \
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

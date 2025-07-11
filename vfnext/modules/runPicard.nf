process runPicard {
  publishDir "${params.outDir}/${sample_id}_results/", mode: "copy", pattern: "{wgs,metrics.alignment*}"
  
  input:
    tuple val(sample_id), path(bams), val(is_paired_end)
    path(ref_fa)

  output:
    tuple val(sample_id), path("wgs"), path("metrics.*")

  script:
  sorted_bam = "${sample_id}.sorted.bam"
  java_cmd = "java -jar /app/picard.jar"
  def count_unpaired = is_paired_end ? "":"--COUNT_UNPAIRED"
  """
  ${java_cmd} CollectWgsMetrics -I ${sorted_bam} \
                                -R ${ref_fa}\
                                -O wgs -CAP 99999 \
                                -Q ${params.base_quality} \
                                -MQ ${params.mapping_quality} \
                                ${count_unpaired}

  ${java_cmd} CollectMultipleMetrics -I ${sorted_bam} \
                                     -R ${ref_fa}\
                                     -O metrics
  """
}

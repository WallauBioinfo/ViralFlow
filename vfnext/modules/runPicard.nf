process runPicard {
  publishDir "${params.outDir}/${sample_id}_results/", mode: "copy", pattern: "{wgs,metrics.alignment*}"
  
  input:
    tuple val(sample_id), path(bams)
    path(ref_fa)

  output:
    tuple val(sample_id), path("wgs"), path("metrics.*")

  script:
  sorted_bam = "${sample_id}.sorted.bam"
  java_cmd = "java -jar /app/picard.jar"
  """
  ${java_cmd} CollectWgsMetrics -I ${sorted_bam} \
                                -R ${ref_fa}\
                                -O wgs -CAP 99999 \
                                -Q ${params.base_quality} \
                                -MQ ${params.mapping_quality}

  ${java_cmd} CollectMultipleMetrics -I ${sorted_bam} \
                                     -R ${ref_fa}\
                                     -O metrics
  """
}

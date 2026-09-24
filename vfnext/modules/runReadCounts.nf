process runReadCounts{
  tag "${meta.id}"
  publishDir { "${params.outDir}/${meta.id}_results/" }, mode: "copy"

  input:
  tuple val(meta), path(bams), val(is_paired_end)
  path(refFa)
  val(depth)
  output:
  tuple val(meta), path("${meta.id}.depth${depth}.fa.bc")

  script:
  sorted_bam = "${bams[0].getSimpleName()}.sorted.bam"
  """
  # RUN READ COUNT
  bam-readcount -d 50000 -q 30 -w 0 -f ${refFa} ${sorted_bam} > ${meta.id}.depth${depth}.fa.bc
  """
}

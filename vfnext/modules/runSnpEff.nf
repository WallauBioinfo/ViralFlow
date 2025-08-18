process runSnpEff{
  tag "${meta.id}"
  errorStrategy 'ignore'
  publishDir "${params.outDir}/${meta.id}_results/", mode: "copy"

  input:
    tuple val(meta), path(bam_files), val(is_paired_end)
    val(genome_code)
    path(refGenomeFasta)
    path(refIndexFiles)
    path(entry_found_file) // here to assure runSnpEff will not start before DB was checked
  
  output:
  tuple val(meta.id), path("*.vcf"), path("${meta.id}_snpEff_summary.html"), path("snpEff_genes.txt")

  script:
  ref_fa = "${refGenomeFasta}"
  sorted_bam = "${bam_files[0].getSimpleName()}.sorted.bam"
  if (params.virus == "custom")
    """
    freebayes -p 1 --reference-quality ${params.mapping_quality},${params.base_quality} \
            -f ${ref_fa} ${sorted_bam} > ${meta.id}.vcf
    snpEff ann -Xmx4g \
            -v ${genome_code} ${meta.id}.vcf > ${meta.id}.ann.vcf
    # add sample id to htmls
    mv snpEff_summary.html ${meta.id}_snpEff_summary.html
    """
  else
    """
    freebayes -p 1 --reference-quality ${params.mapping_quality},${params.base_quality} \
            -f ${ref_fa} ${sorted_bam} > ${meta.id}.vcf
    snpEff ann -Xmx4g \
            ${genome_code} ${meta.id}.vcf > ${meta.id}.ann.vcf
    # add sample id to htmls
    mv snpEff_summary.html ${meta.id}_snpEff_summary.html
    """
}

process runSnpEff{
  publishDir "${params.outDir}/${sample_id}_results/", mode: "copy"

  input:
  val(genome_code)
  path(entry_found_file)
  tuple val(sample_id), path(consensus), path(ivar_txt), path(mut_tsv), path(vcf_file), path(vcf_index)
  
  output:
  tuple val(sample_id), path("*.vcf"), path("${sample_id}_snpEff_summary.html"), path("snpEff_genes.txt")

  script:
  """
  snpEff ann -Xmx4g \
    ${ params.virus == "custom" ? "-v ${genome_code} " : "${genome_code} " } \
    ${vcf_file} > ${sample_id}.ann.vcf 2> snpEff.log
  
  mv snpEff_summary.html ${sample_id}_snpEff_summary.html
  """
}

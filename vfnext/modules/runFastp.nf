
process runFastp{
  publishDir "${params.outDir}/${sample_id}_results/", mode: "copy", pattern: "*.fastp.html"
  label "multithread"

  input:
    tuple val(sample_id), path(reads), val(is_paired_end)

  output:
    // tuple val(prfx), path('*.R{1,2}.fq.gz'), path("${prfx}.fastp.html")
    tuple val(prfx), path('*.*.fq.gz'), path("${prfx}.fastp.html")
    
  script:
      prfx = sample_id
      def dedup = params.dedup ? "--dedup --dup_calc_accuracy ${params.ndedup}":"--dont_eval_duplication"

      """
      if [[ ${is_paired_end}  == true ]]; then
          fastp -i ${reads[0]} -I ${reads[1]} \
            --detect_adapter_for_pe \
            --thread ${params.fastp_threads} \
            -o ${prfx}.R1.fq.gz \
            -O ${prfx}.R2.fq.gz \
            -h ${prfx}.fastp.html \
            -j ${prfx}.fastp.json \
            -l ${params.minLen} -f ${params.trimLen} -t ${params.trimLen} \
            -F ${params.trimLen} -T ${params.trimLen} \
            --cut_front --cut_tail --qualified_quality_phred 20 \
            ${dedup}
      else
          fastp -i ${reads[0]} \
            --thread ${params.fastp_threads} \
            -o ${prfx}.SE.fq.gz \
            -h ${prfx}.fastp.html \
            -j ${prfx}.fastp.json \
            -l ${params.minLen} -f ${params.trimLen} -t ${params.trimLen} \
            -F ${params.trimLen} -T ${params.trimLen} \
            --cut_front --cut_tail --qualified_quality_phred 20 \
            ${dedup}
      fi
      """
}

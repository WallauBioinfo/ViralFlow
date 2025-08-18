process runNextClade {
  publishDir "${params.outDir}/${meta.id}_results/", mode: "copy", pattern: "{*nextclade.csv,*.errors.csv,*.translation.fasta}"
  input:
  tuple val(meta), path(intrahost_tsvs), path(algn_fasta), path(consensus_fa), path(ivar_txt), path(mut_tsv)
  path(ref_fa)

  // temporary solution, no need for ivar_txt and mut_tsv
  output:
  tuple path("*.csv"), path("*.fasta")
  
  shell:
  nxt_dataset = "${workflow.projectDir}/containers/nextclade_dataset/sars-cov-2/"
  '''
  # check if intravariants are present:
  #     if yes, compile all on a single file and use it as input
  #     if no, use the consensus sequence as input

  NUMLINES=$(wc -l < !{meta.id}.depth!{params.depth}.fa.bc.intrahost.short.tsv)

  if [ $NUMLINES -gt 1 ]; then
      cat !{meta.id}.depth!{params.depth}.fa !{meta.id}.depth!{params.depth}.fa.algn.minor.fa  > !{meta.id}.depth!{params.depth}.all.fa
      nextclade run --jobs !{params.nxtclade_jobs} \
                --input-root-seq=!{ref_fa} \
                --input-dataset=!{nxt_dataset} \
                --output-csv=!{meta.id}.depth!{params.depth}.all.fa.nextclade.csv \
                --output-all=./ \
                !{meta.id}.depth!{params.depth}.all.fa
  fi

  if [ $NUMLINES -eq 1 ]; then
      nextclade run --jobs !{params.nxtclade_jobs} \
             --input-root-seq=!{ref_fa} \
             --input-dataset=!{nxt_dataset} \
             --output-csv=!{meta.id}.depth!{params.depth}.fa.nextclade.csv \
             --output-all=./ \
             !{meta.id}.depth!{params.depth}.fa
  fi

  '''
}

process runPangolin {
  publishDir "${params.outDir}/${meta.id}_results/",mode: "copy"
  input:
  tuple val(meta), path(intrahost_tsvs), path(algn_fasta), path(consensus_fa), path(ivar_txt), path(mut_tsv)
  //temporary solution, no need for ivar_txt and mut_tsv

  output:
   path("*.csv")
  shell:
  '''
  NUMLINES=$(wc -l < !{meta.id}.depth!{params.depth}.fa.bc.intrahost.short.tsv)

  if [ $NUMLINES -gt 1 ]; then
      cat !{meta.id}.depth!{params.depth}.fa !{meta.id}.depth!{params.depth}.fa.algn.minor.fa  > !{meta.id}.depth!{params.depth}.all.fa
      pangolin !{meta.id}.depth!{params.depth}.all.fa \
                -t !{params.pangolin_threads} --outfile !{meta.id}.all.fa.pango.out.csv
  fi
  if [ $NUMLINES -eq 1 ]; then
      pangolin !{meta.id}.depth!{params.depth}.fa \
                -t !{params.pangolin_threads} --outfile !{meta.id}.fa.pango.out.csv
  fi

  '''
}

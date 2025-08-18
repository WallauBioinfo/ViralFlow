process runIntraHostScript{
  tag "${meta.id}"
  publishDir "${params.outDir}/${meta.id}_results/", mode: "copy"

  input:
     tuple val(meta), path(fa_bc), path(fa_algn)
     path(ref_gff)

  output:
     tuple val(meta), path("*.tsv"), path("*.fa")

  script:
     """
     python $projectDir/bin/intrahost_scriptv2.py \
            -in ${meta.id}.depth${params.depth}.fa.bc \
            -al ${meta.id}.depth${params.depth}.fa.algn \
            -dp ${params.minDpIntrahost} \
            -gf ${ref_gff}
     """
}

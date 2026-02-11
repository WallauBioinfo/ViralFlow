process alignConsensus2Ref {
    tag "${meta.id}"
    publishDir "${params.outDir}/${meta.id}_results/", mode: "copy"
    label "multithread"

    input:
    tuple val(sample_id), path(consensus_fa), path(ivar_txt), path(mut_tsv), path(vcf_file), path(vcf_index)
    path(ref_fa)
    //temporary solution, we only need the consensus_fa. handle channels better
    output:
    tuple val(meta), path("${meta.id}.depth${params.depth}.fa.algn")
    //path("*.fa.algn")
    script:
    """
    mafft --keeplength --add ${meta.id}.depth${params.depth}.fa \
                       --thread ${params.mafft_threads} \
                       ${ref_fa} \
    > ${meta.id}.depth${params.depth}.fa.algn
    """
}

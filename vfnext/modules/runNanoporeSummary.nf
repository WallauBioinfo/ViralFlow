process runNanoporeSummary {
    label "NP_basecontainer"
    publishDir { "${params.outDir}/${meta.id}_results/" }, mode: 'copy', overwrite: true
    tag "${meta.id}"

    input:
        tuple val(meta),
              path(raw_vcf),
              path(filtered_vcf),
              path(filtered_tbi),
              path(consensus),
              path(low_cov),
              path(coverage)
        path(ref)
        val(clair3_qual)
        val(mapping_quality)
        val(af_threshold)
        val(min_depth)

    output:
        tuple val(meta), path("${meta.id}.nanopore_summary.tsv")

    script:
    // Called by name: Nextflow puts bin/ on the task's PATH wherever the task
    // runs. ${projectDir}/bin names the launch host's copy, which a cloud
    // executor's task never sees; Nextflow uploads bin/ there instead.
    """
    set -euo pipefail

    nanopore_summary.py \
        --raw-vcf ${raw_vcf} \
        --filtered-vcf ${filtered_vcf} \
        --consensus ${consensus} \
        --coverage ${coverage} \
        --reference ${ref} \
        --clair3-qual ${clair3_qual} \
        --mapping-quality ${mapping_quality} \
        --af-threshold ${af_threshold} \
        --min-depth ${min_depth} \
        --output ${meta.id}.nanopore_summary.tsv
    """
}

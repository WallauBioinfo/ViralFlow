process run_nanopore_qc {
    label "NP_basecontainer"
    publishDir { "${workflow.outputDir}/${meta.id}_results/" }, mode: 'copy', overwrite: true
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
        tuple val(meta), path("${meta.id}.nanopore_qc.tsv")

    script:
    """
    set -euo pipefail

    python3 ${projectDir}/bin/nanopore_qc.py \
        --raw-vcf ${raw_vcf} \
        --filtered-vcf ${filtered_vcf} \
        --consensus ${consensus} \
        --coverage ${coverage} \
        --reference ${ref} \
        --clair3-qual ${clair3_qual} \
        --mapping-quality ${mapping_quality} \
        --af-threshold ${af_threshold} \
        --min-depth ${min_depth} \
        --output ${meta.id}.nanopore_qc.tsv
    """
}

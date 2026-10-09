
process runBcftools {
    label "NP_basecontainer"
    // Define the process parameters
    publishDir { "${params.outDir}/${meta.id}_results/" }, mode: 'copy', overwrite: true
    tag "${meta.id}"

    input:
        tuple val(meta), path(vcf)
        path(ref)
        val(af_threshold) // 0.51

    output:
        tuple val(meta),
              path("${meta.id}.filtered.vcf.gz"),
              path("${meta.id}.filtered.vcf.gz.tbi")

    script:
    // FILTER="PASS" as well as the AF cutoff: Clair3's --qual (clair3_qual)
    // does not drop a call below it, it labels it LowQual and keeps it. Without
    // this, LowQual calls reached the consensus and clair3_qual changed only
    // the summary. Clair3's RefCall rows, printed only with --print_ref_calls,
    // are not PASS either.
    """
    set -euo pipefail

    bcftools norm -m - -f ${ref} ${vcf} -Ou \
        | bcftools filter -i 'FILTER="PASS" && FORMAT/AF >= ${af_threshold}' -Oz \
            -o ${meta.id}.filtered.vcf.gz

    bcftools index --tbi ${meta.id}.filtered.vcf.gz
    """
}

process runBcftoolsConsensus {
    label "NP_basecontainer"
    // Define the process parameters
    publishDir { "${params.outDir}/${meta.id}_results/" }, mode: 'copy', overwrite: true
    tag "${meta.id}"

    input:
        tuple val(meta), path(vcf), path(tbi), path(bam), path(bai)
        path(ref)
        val(min_depth) // 3

    output:
        tuple val(meta), path("${meta.id}.consensus.fa"), path("${meta.id}.low_cov.bed"), path("${meta.id}.cov.bed")

    script:
    """
    set -euo pipefail

    # create a bed file with low coverage regions
    # this is used to mask low coverage regions in the consensus sequence
    #
    # -aa, not -a: a single -a reports zero-depth positions only on contigs
    # that have at least one read, and says nothing at all about a contig no
    # read reached. Such a contig was therefore never masked, so a sample with
    # no aligned reads - a negative control, a failed barcode - published the
    # reference itself as its consensus, reported 100% callable.
    samtools depth -J -aa ${bam} > ${meta.id}.cov.bed
    awk '\$3 <= int(${min_depth}) {print \$1 "\t" \$2-1 "\t" \$2}' ${meta.id}.cov.bed > ${meta.id}.low_cov.bed

    # call consensus sequence and rename it
    bcftools consensus -f ${ref} --mask ${meta.id}.low_cov.bed ${vcf} > ${meta.id}.consensus.fa
    sed -i -e 's/>.*/>${meta.id}/' ${meta.id}.consensus.fa
    """
}

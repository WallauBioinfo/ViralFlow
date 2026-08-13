nextflow.enable.dsl = 2

include {
    run_bcftools
    run_bcftools_consensus
} from '../../modules/runBcftools.nf'
include { run_nanopore_qc } from '../../modules/runNanoporeQc.nf'

process prepare_fixture_bam {
    label "NP_basecontainer"
    tag "${meta.id}"

    input:
        tuple val(meta), path(sam)

    output:
        tuple val(meta),
              path("${meta.id}.sorted.bam"),
              path("${meta.id}.sorted.bam.bai")

    script:
    """
    set -euo pipefail

    samtools view -bS ${sam} \
        | samtools sort -o ${meta.id}.sorted.bam

    samtools index ${meta.id}.sorted.bam
    """
}

workflow BCFTOOLS_FIXTURE {
    take:
        vcf_ch
        ref
        sam_ch

    main:
        prepare_fixture_bam(sam_ch)
        run_bcftools(vcf_ch, ref, 0.51)

        prepare_fixture_bam.out
            .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
            .set { keyed_bams_ch }

        run_bcftools.out
            .map { meta, vcf, tbi -> tuple(meta.id, vcf, tbi) }
            .set { keyed_vcfs_ch }

        keyed_bams_ch
            .join(keyed_vcfs_ch)
            .map { _id, meta, bam, bai, vcf, tbi ->
                tuple(meta, vcf, tbi, bam, bai)
            }
            .set { consensus_input_ch }

        run_bcftools_consensus(consensus_input_ch, ref, params.fixture_min_depth)

        vcf_ch
            .map { meta, vcf -> tuple(meta.id, vcf) }
            .set { keyed_raw_vcfs_ch }

        run_bcftools_consensus.out
            .map { meta, consensus, low_cov, coverage ->
                tuple(meta.id, meta, consensus, low_cov, coverage)
            }
            .set { keyed_consensus_ch }

        keyed_raw_vcfs_ch
            .join(keyed_vcfs_ch)
            .join(keyed_consensus_ch)
            .map { _id, raw_vcf, filtered_vcf, filtered_tbi, meta, consensus, low_cov, coverage ->
                tuple(meta, raw_vcf, filtered_vcf, filtered_tbi, consensus, low_cov, coverage)
            }
            .set { qc_input_ch }

        run_nanopore_qc(
            qc_input_ch,
            ref,
            10,
            30,
            0.51,
            params.fixture_min_depth
        )

    emit:
        filtered = run_bcftools.out
        consensus = run_bcftools_consensus.out
        qc = run_nanopore_qc.out
}

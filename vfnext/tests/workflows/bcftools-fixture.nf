nextflow.enable.dsl = 2

include {
    runBcftools
    runBcftoolsConsensus
} from '../../modules/runBcftools.nf'
include { runNanoporeSummary } from '../../modules/runNanoporeSummary.nf'

process prepareFixtureBam {
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
        vcfCh
        ref
        samCh

    main:
        prepareFixtureBam(samCh)
        runBcftools(vcfCh, ref, 0.51)

        prepareFixtureBam.out
            .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
            .set { keyedBamsCh }

        runBcftools.out
            .map { meta, vcf, tbi -> tuple(meta.id, vcf, tbi) }
            .set { keyedVcfsCh }

        keyedBamsCh
            .join(keyedVcfsCh)
            .map { _id, meta, bam, bai, vcf, tbi ->
                tuple(meta, vcf, tbi, bam, bai)
            }
            .set { consensusInputCh }

        runBcftoolsConsensus(consensusInputCh, ref, params.fixture_min_depth)

        vcfCh
            .map { meta, vcf -> tuple(meta.id, vcf) }
            .set { keyedRawVcfsCh }

        runBcftoolsConsensus.out
            .map { meta, consensus, low_cov, coverage ->
                tuple(meta.id, meta, consensus, low_cov, coverage)
            }
            .set { keyedConsensusCh }

        keyedRawVcfsCh
            .join(keyedVcfsCh)
            .join(keyedConsensusCh)
            .map { _id, raw_vcf, filtered_vcf, filtered_tbi, meta, consensus, low_cov, coverage ->
                tuple(meta, raw_vcf, filtered_vcf, filtered_tbi, consensus, low_cov, coverage)
            }
            .set { qcInputCh }

        runNanoporeSummary(
            qcInputCh,
            ref,
            10,
            30,
            0.51,
            params.fixture_min_depth
        )

    emit:
        filtered = runBcftools.out
        consensus = runBcftoolsConsensus.out
        summary = runNanoporeSummary.out
}

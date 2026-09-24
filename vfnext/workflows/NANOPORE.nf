//enable dsl 2
nextflow.enable.dsl = 2
include {runPorechop} from '../modules/runPorechop.nf'
include {runMinimap2} from '../modules/runMinimap2.nf'
include {runFaidx} from '../modules/runFaidx.nf'
include {runAmpliconClip} from '../modules/runAmpliconClip.nf'
include {runBamUtils} from '../modules/runBamUtils.nf'
include {runClair3} from '../modules/runClair3.nf'
include {runBcftools; runBcftoolsConsensus} from '../modules/runBcftools.nf'
include {runNanoporeSummary} from '../modules/runNanoporeSummary.nf'
include {normalizeTrimLen} from '../modules/param_helpers.nf'

workflow NANOPORE {
    take:
        readsCh // tuple (meta, fastq)
        ref // path to reference genome

    main:

    // remove adapters (porechop)
    runPorechop(readsCh)
    readsCh = runPorechop.out

    // do alignment (minimap2)
    runMinimap2(readsCh, ref)
    bamsCh = runMinimap2.out // tuple (meta, sorted_bam, bai)

    // optional primer clipping (samtools ampliconclip)
    //
    // Off unless --primersBED is supplied, which matches ILLUMINA. When it is
    // supplied the clipped BAM replaces the raw alignment for everything
    // downstream, so variant calling, the depth used for masking and the
    // consensus all see primer-free reads. Leaving it out of only some of those
    // would make the consensus disagree with the variants it was built from.
    if (params.primersBED) {
        runAmpliconClip(bamsCh, file(params.primersBED))
        bamsCh = runAmpliconClip.out.bams
    }

    // optional read-end trimming (bamUtil trimBam)
    //
    // Off while trimLen is 0, its default. ILLUMINA applies the same parameter
    // in fastp, before alignment; there is no equivalent FASTQ step here, so
    // NANOPORE trims the aligned BAM instead. Runs after primer clipping so the
    // trim removes bases beyond the primers rather than eating into the region
    // ampliconclip is about to look for.
    def trim_len = normalizeTrimLen(params.trimLen)
    if (trim_len > 0) {
        runBamUtils(bamsCh, trim_len)
        bamsCh = runBamUtils.out.bams
    }

    // index the reference once, not inside every Clair3 task
    //
    // There is a single reference per run, so runFaidx emits one .fai and
    // Nextflow broadcasts it across the per-sample BAM channel. No .first() or
    // .collect() is needed, and adding one would be actively wrong if the
    // reference ever became per-sample: it would silently pin every sample to
    // the first index instead of failing.
    runFaidx(ref)

    // do variant calling (clair3)
    runClair3(bamsCh, ref, runFaidx.out.fai, params.clair3_chunk_size, params.clair3_qual, params.mapping_quality, params.clair3_model)

    // normlalize indes and filter variants (bcftools)
    runBcftools(runClair3.out, ref, params.af_threshold)

    bamsCh
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

    // call consensus sequence (bcftools consensus)
    runBcftoolsConsensus(consensusInputCh, ref, params.np_min_depth)

    runClair3.out
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
        .set { nanoporeSummaryInputCh }

    runNanoporeSummary(
        nanoporeSummaryInputCh,
        ref,
        params.clair3_qual,
        params.mapping_quality,
        params.af_threshold,
        params.np_min_depth
    )

    bamsCh
        .map { meta, bam, bai -> tuple(meta, bam, bai, false) }
        .set { plotBamsCh }

    emit:
        bamsCh = plotBamsCh // tuple (meta, sorted_bam, bai, is_paired_end)
        rawVcfsCh = runClair3.out
        filteredVcfsCh = runBcftools.out
        consensusCh = runBcftoolsConsensus.out
        summaryCh = runNanoporeSummary.out
}

//enable dsl 2
nextflow.enable.dsl = 2
include {run_porechop} from '../modules/runPorechop.nf'
include {run_minimap2} from '../modules/runMinimap2.nf'
include {run_amplicon_clip} from '../modules/runAmpliconClip.nf'
include {run_bam_utils} from '../modules/runBamUtils.nf'
include {run_clair3} from '../modules/runClair3.nf'
include {run_bcftools; run_bcftools_consensus} from '../modules/runBcftools.nf'
include {run_nanopore_qc} from '../modules/runNanoporeQc.nf'

// nextflow.config declares trimLen as an Integer, but a value given on the
// command line arrives as a String, and the wrapper forwards params-file
// entries as command line arguments too. Comparing the two aborts the run
// before any task is submitted ("Cannot compare java.lang.String with value
// '30' and java.lang.Integer with value '0'"), so normalize once, here rather
// than in step0 validation: the integration tests call NANOPORE directly,
// without processInputs.
def normalizeTrimLen(value) {
    def raw = value == null ? '0' : value.toString().trim()
    if (!(raw ==~ /\d+/)) {
        error "--trimLen must be a non-negative integer, got '${value}'"
    }
    return raw as Integer
}

workflow NANOPORE {
    take:
        reads_ch // tuple (meta, fastq)
        ref // path to reference genome

    main:

    // remove adapters (porechop)
    run_porechop(reads_ch)
    reads_ch = run_porechop.out

    // do alignment (minimap2)
    run_minimap2(reads_ch, ref)
    bams_ch = run_minimap2.out // tuple (meta, sorted_bam, bai)

    // optional primer clipping (samtools ampliconclip)
    //
    // Off unless --primersBED is supplied, which matches ILLUMINA. When it is
    // supplied the clipped BAM replaces the raw alignment for everything
    // downstream, so variant calling, the depth used for masking and the
    // consensus all see primer-free reads. Leaving it out of only some of those
    // would make the consensus disagree with the variants it was built from.
    if (params.primersBED) {
        run_amplicon_clip(bams_ch, file(params.primersBED))
        bams_ch = run_amplicon_clip.out.bams
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
        run_bam_utils(bams_ch, trim_len)
        bams_ch = run_bam_utils.out.bams
    }

    // do variant calling (clair3)
    run_clair3(bams_ch, ref, params.clair3_chunk_size, params.clair3_qual, params.mapping_quality, params.clair3_model)

    // normlalize indes and filter variants (bcftools)
    run_bcftools(run_clair3.out, ref, params.af_threshold)

    bams_ch
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

    // call consensus sequence (bcftools consensus)
    run_bcftools_consensus(consensus_input_ch, ref, params.np_min_depth)

    run_clair3.out
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
        .set { nanopore_qc_input_ch }

    run_nanopore_qc(
        nanopore_qc_input_ch,
        ref,
        params.clair3_qual,
        params.mapping_quality,
        params.af_threshold,
        params.np_min_depth
    )

    bams_ch
        .map { meta, bam, bai -> tuple(meta, bam, bai, false) }
        .set { plot_bams_ch }

    emit:
        bams_ch = plot_bams_ch // tuple (meta, sorted_bam, bai, is_paired_end)
        raw_vcfs_ch = run_clair3.out
        filtered_vcfs_ch = run_bcftools.out
        consensus_ch = run_bcftools_consensus.out
        qc_ch = run_nanopore_qc.out
}

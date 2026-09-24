
include { runAmpliconClipping } from '../modules/ampliconclip.nf'
include { align2ref } from '../modules/align2ref.nf'
include { runIvar } from '../modules/runIvar.nf'
include { checkSnpEffDB } from '../modules/checkSnpEffDB.nf'
include { runSnpEff } from '../modules/runSnpEff.nf'
include { runVfReport } from '../modules/runVfReport.nf'
include { alignConsensus2Ref } from '../modules/alignConsensus2Ref.nf'
include { snpPlot } from '../modules/generatePlots.nf'
include { runPicard } from '../modules/runPicard.nf'
include { fixWGS } from '../modules/fixWGS.nf'
include { runIntraHostScript } from '../modules/runIntraHostScript.nf'
include { runPangolin } from '../modules/runPangolin.nf'
include { runNextClade } from '../modules/runNextclade.nf'
include { compileOutputs; compileOutputs as compileOutputs_SC2} from '../modules/compileOutput.nf'
include { runReadCounts } from '../modules/runReadCounts.nf'
include { genFaIdx } from '../modules/genFaIdx.nf'
include { indexReferenceBWA } from '../modules/bwaIndex.nf'
include { runFastp } from '../modules/runFastp.nf'

workflow  ILLUMINA {
    take:
        readsCh
        refFa
        refGff
        refGcode
    main:

    // STEP 1 ------------------------------------------------------------------
    // run indexing, open the bwa index output channel
    indexReferenceBWA(refFa)
    indexReferenceBWA.out.set { bwaidxOutputCh }


    // run fastp
    readsCh2 = readsCh.map { meta, files ->
    tuple(meta, files, meta.is_paired_end)  // Add is paired end to the tuple
    }

    runFastp(readsCh2)

    // collect htmls for vf reports
    runFastp.out //tuple (meta, [fq.gz file(s)], fastp_html)
      | map { output -> output[2] }
      | set {fastpHtmlCh}
    allFastpHtmlCh = fastpHtmlCh.collect()

    // collect output reads
    runFastp.out // tuple (meta, [fq.gz file(s)], fastp_html)
      | map { output -> tuple(output[0], output[1]) } //tuple (meta, [fq.gz file(s)])
      | set {fastpFqgzCh}

    fastpFqgzCh = fastpFqgzCh.map { meta, files ->
    def is_paired_end = (files.size() == 2) // Check if it's paired-end
    def new_meta = meta.plus([is_paired_end: is_paired_end]) // Add is_paired_end to the metadata
    tuple(new_meta, files, is_paired_end)  // Add is paired end to the tuple
    }

    // generate fa index
    genFaIdx(refFa)

    // align 2 reference -----------------------------------------------------------
    align2refInCh = fastpFqgzCh.combine(bwaidxOutputCh) // tuple(meta, reads, is_paired_end, fasta_amb, fasta_ann, fasta_bwt, fasta_pac, fasta_sa)

    align2ref(align2refInCh, refFa)
    // Conditionally run ampliconclip for primer trimming if BED file is provided
    if (params.primersBED != null) {
      runAmpliconClipping(align2ref.out.regular_output,
                          params.primersBED,
                          params.minLen)
      bamOutputCh = runAmpliconClipping.out.regular_output
    } else {
      bamOutputCh = align2ref.out.regular_output
    }

    // use bam output for downstream processing
    bamOutputCh
      | map { meta, bam, _bai, is_pe -> tuple(meta, bam, is_pe) }
      | set { bamOutCh }

    // remove bam files which are too small (necessary for Picard)
    bamOutCh
      | filter { output ->
        def bam = output[1]
        // filter if unix paths
        ((bam.getClass() == sun.nio.fs.UnixPath) && (bam.size() >= params.minBamSize ))
        ||
        ((bam.getClass() == java.util.ArrayList) && (bam[0].size() >= params.minBamSize ) && (bam[1].size() >= params.minBamSize))
      }
      | set {bamOutFilteredCh}

    // raise warning in case anyfile is excluded
    bamOutCh
      | filter { output ->
        def bam = output[1]
        // filter if unix paths
        ((bam.getClass() == sun.nio.fs.UnixPath) && (bam.size() <= params.minBamSize ))
        ||
        // filter if is a list with two bam file paths
        ((bam.getClass() == java.util.ArrayList) && (bam[0].size() <= params.minBamSize ) && (bam[1].size() <= params.minBamSize))
      }
      | view { output -> log.warn("Excluding ${output[0].id} bam files as input for Picard due to small size (< ${params.minBamSize} bytes)") }


    // -----------------------------------------------------------------------------
    // Call consensus
    // ivar
    runIvar(bamOutCh, refFa)
    runIvar.out.set { runIvarOutCh }

    // get VCFs
    if ((params.runSnpEff==true)) {
      // check if genome code is on SnpEff database
      checkSnpEffDB(refGcode)
      // runSnpEffDB
      runSnpEff(refGcode,
		checkSnpEffDB.out,
		runIvarOutCh)

	runSnpEff.out
		| map { output -> output[2] }
		| set { snpEff_html }

        allSnpEffHtmlCh = snpEff_html.collect()
        runVfReport(allFastpHtmlCh, allSnpEffHtmlCh)
    }


    // align consensus to ref
    alignConsensus2Ref(runIvarOutCh, refFa)
    alignConsensus2Ref.out.set {alignConOutCh}

    // Rendering snp plot
    snpPlot(alignConOutCh)

    // Assembly Metrics
    runPicard(bamOutFilteredCh, refFa)
    runPicard.out.set {runPicardOutCh}
    fixWGSInCh = runPicardOutCh.join(runIvarOutCh)
    fixWGS(fixWGSInCh)
    // readcounts
    runReadCounts(bamOutCh, refFa, params.depth)
    runReadCounts.out.set {runReadCountsOutCh}
    // run intrahost
    intraHostInCh = alignConOutCh.join(runReadCountsOutCh)
    runIntraHostScript(intraHostInCh, refGff)
    runIntraHostScript.out.set {runIntraHostScriptOutCh}

    // run Variant Naming (Pangolin and Nextclade)
    // runIvar also emits the compressed VCF and its index.  The lineage tools
    // only consume the consensus-related fields, so keep their input tuple at
    // the six elements declared by runPangolin and runNextClade.
    runVariantNamingInCh = runIntraHostScriptOutCh
        .join(runIvarOutCh)
        .map { meta, intrahost_tsvs, algn_fasta, consensus_fa, ivar_txt, mut_tsv, _vcf_file, _vcf_index ->
            tuple(meta, intrahost_tsvs, algn_fasta, consensus_fa, ivar_txt, mut_tsv)
        }

    if (params.virus=="sars-cov2"){
        runPangolin(runVariantNamingInCh)
        runNextClade(runVariantNamingInCh, refFa)
        // Pangolin is the last ones to run, so will use it as a trigger to
        // the output compilation/
        // for the final version, need to find a better way. Maybe split and set
        // as individual post analysis workflow
        finalTrigger = runPangolin.out.concat(fixWGS.out).collect()
        compileOutputs_SC2(finalTrigger, params.virus)
    }

    if (params.virus=="custom"){
        // GAMBIARRA ALERT
        finalTrigger = runIntraHostScript.out.concat(fixWGS.out).collect()
        compileOutputs(finalTrigger, params.virus)
  }

  emit:
    bamsCh = bamOutputCh // meta, sorted_bam, bai, is_paired_end
}

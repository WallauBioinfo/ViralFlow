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
        reads_ch
        ref_fa
        ref_gff
        ref_gcode
    main:

  // STEP 1 ------------------------------------------------------------------
  // run indexing, open the bwa index output channel
  indexReferenceBWA(ref_fa)
  indexReferenceBWA.out.set { bwaidx_Output_ch }


  // run fastp
  reads_ch_2 = reads_ch.map { meta, files ->
  tuple(meta, files, meta.is_paired_end)  // Add is paired end to the tuple
  }

  runFastp(reads_ch_2)

  // collect htmls for vf reports
  runFastp.out //tuple (meta, [fq.gz file(s)], fastp_html)
    | map{ it[2]} 
    | set {fastp_html_ch}
  all_fastp_html_ch = fastp_html_ch.collect()
  
  // collect output reads
  runFastp.out // tuple (meta, [fq.gz file(s)], fastp_html)
    | map {tuple(it[0],it[1])} //tuple (meta, [fq.gz file(s)])
    | set {fastp_fqgz_ch}

  fastp_fqgz_ch = fastp_fqgz_ch.map { meta, files ->
   def is_paired_end = (files.size() == 2) // Check if it's paired-end
   def new_meta = meta.plus([is_paired_end: is_paired_end]) // Add is_paired_end to the metadata
   tuple(new_meta, files, is_paired_end)  // Add is paired end to the tuple
  }


    // generate fa index
    genFaIdx(ref_fa)
    genFaIdx.out.set {faIdx_ch}

    // align 2 reference -----------------------------------------------------------
    align2ref_In_ch = fastp_fqgz_ch.combine(bwaidx_Output_ch) // tuple(meta, reads, is_paired_end, fasta_amb, fasta_ann, fasta_bwt, fasta_pac, fasta_sa)
   
    align2ref(align2ref_In_ch, ref_fa)
    // remove bai file (not used downstream, but usefull as a pipeline output)
    align2ref.out.regular_output // tuple (meta, bam_file, bai_file, is paired_end)
    | map { tuple(it[0], it[1], it[3]) } // tuple(meta, bam_file, is_paired_end)
    | set { align2ref_Out_ch }

    // remove bam files which are too small (necessary for Picard)
    align2ref_Out_ch
    | filter{
        // filter if unix paths
        ((it[1].getClass() == sun.nio.fs.UnixPath) && (it[1].size() >= params.minBamSize )) 
        ||
        ((it[1].getClass() == java.util.ArrayList) && (it[1][0].size() >= params.minBamSize ) && (it[1][1].size() >= params.minBamSize))
    }
    | set {align2ref_Out_filtered_ch}

    // raise warning in case anyfile is excluded
    align2ref_Out_ch
    | filter{
      // filter if unix paths
      (it[1].getClass() == sun.nio.fs.UnixPath) && (it[1].size() <= params.minBamSize )
    }
    | filter {
      // filter if is a list with two bam file paths
      (it[1].getClass() == java.util.ArrayList) && (it[1][0].size() <= params.minBamSize ) && (it[1][1].size() <= params.minBamSize)
    }
    | view{log.warn("Excluding ${it[0]} bam files as input for Picard due to small size (< ${params.minBamSize} bytes)")}

    // -----------------------------------------------------------------------------
    // Call consensus  
    // ivar
    runIvar(align2ref_Out_ch, ref_fa)
    runIvar.out.set { runIvar_Out_ch }

    // get VCFs
    if ((params.runSnpEff==true)) {
        // check if genome code is on SnpEff database
        checkSnpEffDB(ref_gcode)
        // runSnpEffDB
        runSnpEff(align2ref_Out_ch,
              ref_gcode,
              ref_fa,
              faIdx_ch,
              checkSnpEffDB.out)
    
        runSnpEff.out
        | map {it[2]}
        | set { snpEff_html }
    
        all_snpEff_html_ch = snpEff_html.collect()
        runVfReport(all_fastp_html_ch, all_snpEff_html_ch)
    }


    // align consensus to ref
    alignConsensus2Ref(runIvar_Out_ch, ref_fa)
    alignConsensus2Ref.out.set {alignCon_Out_ch}

    // Rendering snp plot
    snpPlot(alignCon_Out_ch)

    // Assembly Metrics
    runPicard(align2ref_Out_filtered_ch, ref_fa)
    runPicard.out.set {runPicard_Out_ch}
    fixWGS_In_ch = runPicard_Out_ch.join(runIvar_Out_ch)
    fixWGS(fixWGS_In_ch)
    // readcounts
    runReadCounts(align2ref_Out_ch, ref_fa, params.depth)
    runReadCounts.out.set {runReadCounts_Out_ch}
    // run intrahost
    intraHost_In_ch = alignCon_Out_ch.join(runReadCounts_Out_ch)
    runIntraHostScript(intraHost_In_ch, ref_gff)
    runIntraHostScript.out.set {runIntraHostScript_Out_ch}

    // run Variant Naming (Pangolin and Nextclade)
    runVariantNaming_In_ch = runIntraHostScript_Out_ch.join(runIvar_Out_ch)

    if (params.virus=="sars-cov2"){
        runPangolin(runVariantNaming_In_ch)
        runNextClade(runVariantNaming_In_ch, ref_fa)
        // Pangolin is the last ones to run, so will use it as a trigger to
        // the output compilation/
        // for the final version, need to find a better way. Maybe split and set
        // as individual post analysis workflow
        final_trigger = runPangolin.out.concat(fixWGS.out).collect()
        compileOutputs_SC2(final_trigger, params.virus)
    }

    if (params.virus=="custom"){
        // GAMBIARRA ALERT
        final_trigger = runIntraHostScript.out.concat(fixWGS.out).collect()
        compileOutputs(final_trigger, params.virus)
  }

  emit:
    bams_ch = align2ref.out.regular_output // sample_id, sorted_bams, .bais, is_paired_end, 
}

#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

<<<<<<< HEAD
=======
// import modules
include { indexReferenceBWA } from './modules/bwaIndex.nf'
include { runFastp } from './modules/runFastp.nf'
include { align2ref } from './modules/align2ref.nf'
include { runAmpliconClipping } from './modules/ampliconclip.nf'
include { coveragePlot } from './modules/generatePlots.nf'
include { snpPlot } from './modules/generatePlots.nf'
include { runIvar } from './modules/runIvar.nf'
include { runReadCounts } from './modules/runReadCounts.nf'
include { alignConsensus2Ref } from './modules/alignConsensus2Ref.nf'
include { runIntraHostScript } from './modules/runIntraHostScript.nf'
include { runPangolin } from './modules/runPangolin.nf'
include { runNextClade } from './modules/runNextclade.nf'
include { runPicard } from './modules/runPicard.nf'
include { fixWGS } from './modules/fixWGS.nf'
include { compileOutputs } from './modules/compileOutput.nf'
include { compileOutputs as compileOutputs_SC2} from './modules/compileOutput.nf'
include { runSnpEff } from './modules/runSnpEff.nf'
include { genFaIdx } from './modules/genFaIdx.nf'
include { getMappedReads } from './modules/getMappedReads.nf'
include { getUnmappedReads } from './modules/getUnmappedReads.nf'
include { checkSnpEffDB } from './modules/checkSnpEffDB.nf'
>>>>>>> develop
// import sub workflows
include { processInputs } from './workflows/step0-input-handling.nf'
include {ILLUMINA} from './workflows/ILLUMINA.nf'
include {NANOPORE} from './workflows/NANOPORE.nf'
include {GENPLOTS} from './workflows/GENPLOTS.nf'

// The code for the inital log info is based on the one found at FASTQC PIPELINE
// https://github.com/angelovangel/nxf-fastqc/blob/master/main.nf


//  The default workflow
workflow {

/*
* ANSI escape codes to color output messages, get date to use in results folder name
*/
def ANSI_GREEN = "\033[1;32m"
def ANSI_RED = "\033[1;31m"
def ANSI_RESET = "\033[0m"

log.info """
  ===========================================
<<<<<<< HEAD
  VFNEXT ${workflow.manifest.version}
=======
  VFNEXT v1.3.5
>>>>>>> develop
  parameters:
  -------------------------------------------
  --inDir            : ${params.inDir}
  --outDir           : ${params.outDir}
  --virus            : ${params.virus}
  --refGenomeCode   *: ${params.refGenomeCode}
  --referenceGenome *: ${params.referenceGenome}
  --referenceGFF    *: ${params.referenceGFF}
  --primersBED       : ${params.primersBED}
  --minLen           : ${params.minLen}
  --depth            : ${params.depth}
  --minDpIntrahost   : ${params.minDpIntrahost}
  --trimLen          : ${params.trimLen}
  --databaseDir      : ${params.databaseDir}
  --runSnpEff        : ${params.runSnpEff}
  --writeMappedReads : ${params.writeMappedReads}
  --nextflowSimCalls : ${params.nextflowSimCalls}
  --fastp_threads    : ${params.fastp_threads}
  --dedup            : ${params.dedup}
  --ndedup           : ${params.ndedup}
  --bwa_threads      : ${params.bwa_threads}
  --mafft_threads    : ${params.mafft_threads}
  --mapping_quality  : ${params.mapping_quality}
  --base_quality     : ${params.base_quality}
  --minBamSize       : ${params.minBamSize}
         
        
  * Only required for "custom" virus
  Runtime data:
  -------------------------------------------
  Running with profile:   ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
  Running as user:        ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
  Launch dir:             ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
  Base dir:               ${ANSI_GREEN}${baseDir}${ANSI_RESET}
  ------------------------------------------
  """.stripIndent()

  // open input channels
  processInputs()
  reads_ch = processInputs.out.reads_ch
  ref_gff = processInputs.out.ref_gff
  ref_fa = processInputs.out.ref_fa
  ref_gcode = processInputs.out.ref_gcode

<<<<<<< HEAD
  if (params.mode == "ILLUMINA"){
    ILLUMINA(reads_ch, ref_fa,ref_gff,ref_gcode)
    def bams_ch = ILLUMINA.out.bams_ch.map{meta, bams, bais, _is_paired_end -> tuple(meta, bams,bais)}
    GENPLOTS(bams_ch)
=======
  // STEP 1 ------------------------------------------------------------------
  // run indexing, open the bwa index output channel
  indexReferenceBWA(ref_fa)
  indexReferenceBWA.out.set { bwaidx_Output_ch }
  genFaIdx(ref_fa)
  genFaIdx.out.set {faIdx_ch}

  // run fastp
  reads_ch = reads_ch.map { sample_id, files ->
  def is_paired_end = (files.size() == 2) // Check if it's paired-end
  tuple(sample_id, files, is_paired_end)  // Add is paired end to the tuple
  }
  runFastp(reads_ch)


  // collect htmls for vf reports
  runFastp.out //tuple (filename_prefix, [fq.gz file(s)], fastp_html)
    | map{ it -> it[2]} 
    | set {fastp_html_ch}
  all_fastp_html_ch = fastp_html_ch.collect()
  
  // collect output reads
  runFastp.out // tuple (filename_prefix, [fq.gz file(s)], fastp_html)
    | map {it -> tuple(it[0],it[1])} //tuple (filename_prefix, [fq.gz file(s)])
    | set {fastp_fqgz_ch}

  fastp_fqgz_ch = fastp_fqgz_ch.map { sample_id, files ->
   def is_paired_end = (files.size() == 2) // Check if it's paired-end
   tuple(sample_id, files, is_paired_end)  // Add is paired end to the tuple
  }
  
  //align 2 reference
  align2ref_In_ch = fastp_fqgz_ch.combine(bwaidx_Output_ch)
   
  align2ref(align2ref_In_ch, ref_fa)

  // Conditionally run ampliconclip for primer trimming if BED file is provided
  if (params.primersBED != null) {
    runAmpliconClipping(align2ref.out.regular_output,
                        params.primersBED, 
                        params.minLen)
    bam_output_ch = runAmpliconClipping.out.regular_output
  } else {
    bam_output_ch = align2ref.out.regular_output
  }

  // use bam output for downstream processing
  bam_output_ch // tuple (sample_id, bam_file, bai_file, is_paired_end)
    | map { it -> tuple(it[0], it[1], it[3]) } // tuple(sample_id, bam_file, is_paired_end)
    | set { bam_Out_ch }

  // remove bam files which are too small (necessary for Picard)
  bam_Out_ch
    | filter(it -> 
      // filter if unix paths
      ((it[1].getClass() == sun.nio.fs.UnixPath) && (it[1].size() >= params.minBamSize )) 
      ||
      ((it[1].getClass() == java.util.ArrayList) && (it[1][0].size() >= params.minBamSize ) && (it[1][1].size() >= params.minBamSize))
    )
    | set {bam_Out_filtered_ch}

  // raise warning in case anyfile is excluded
  bam_Out_ch
    | filter(it -> 
      // filter if unix paths
      (it[1].getClass() == sun.nio.fs.UnixPath) && (it[1].size() <= params.minBamSize )
    )
    | filter (it ->
      // filter if is a list with two bam file paths
      (it[1].getClass() == java.util.ArrayList) && (it[1][0].size() <= params.minBamSize ) && (it[1][1].size() <= params.minBamSize)
     )
    | view(it -> log.warn("Excluding ${it[0]} bam files as input for Picard due to small size (< ${params.minBamSize} bytes)"))

  //Rendering the depth coverage plot
  coveragePlot(bam_output_ch)
  // Check if there are mapped reads
  coveragePlot_out_ch = coveragePlot.out.result
  coveragePlot_out_ch
    | view(it -> log.warn("${it.text}"))

  if ((params.writeMappedReads == true)){
    // write mapped reads
    getMappedReads(bam_Out_ch)
  
    // write unmappped reads
    getUnmappedReads(bam_Out_ch)
  }
  
  // ivar
  runIvar(bam_Out_ch, ref_fa)
  runIvar.out.set { runIvar_Out_ch }

  // readcounts
  runReadCounts(bam_Out_ch, ref_fa)
  runReadCounts.out.set {runReadCounts_Out_ch}

  // get VCFs
  if ((params.runSnpEff==true)) {
    // check if genome code is on SnpEff database
    checkSnpEffDB(ref_gcode)
    // runSnpEffDB
    runSnpEff(ref_gcode,
              checkSnpEffDB.out,
              runIvar_Out_ch)
    
    runSnpEff.out
      | map {it -> it[2]}
      | set { snpEff_html }
    
    all_snpEff_html_ch = snpEff_html.collect()
    runVfReport(all_fastp_html_ch, all_snpEff_html_ch)
  }
  
  //align consensus to ref
  alignConsensus2Ref(runIvar_Out_ch, ref_fa)
  alignConsensus2Ref.out.set {alignCon_Out_ch}

  //Rendering snp plot
  snpPlot(alignCon_Out_ch)

  // Assembly Metrics
  runPicard(bam_Out_filtered_ch, ref_fa)
  runPicard.out.set {runPicard_Out_ch}
  fixWGS_In_ch = runPicard_Out_ch.join(runIvar_Out_ch)
  fixWGS(fixWGS_In_ch)

  //run intrahost
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
>>>>>>> develop
  }

  if (params.mode == "NANOPORE"){
    NANOPORE(reads_ch, ref_fa)
    def bams_ch = NANOPORE.out.bams_ch.map{meta, bam -> 
        // Remove once NANOPORE emits (meta,bam,bai)
        assert bam instanceof List && bam.size()==2 : "NANOPORE bams_ch must be meta, [bam,bai]"
        tuple(meta, bam[0], bam[1])}
    GENPLOTS(bams_ch)
  }

// -------------- Check if everything went okay -------------------------------
workflow.onComplete = {
  if (workflow.success) {
    log.info """
      ===========================================
      ${ANSI_GREEN}Finished in ${workflow.duration}
      """.stripIndent()
  } else {
    log.info """
      ===========================================
      ${ANSI_RED}Finished with errors!${ANSI_RESET}
      """.stripIndent()
  }
}

}


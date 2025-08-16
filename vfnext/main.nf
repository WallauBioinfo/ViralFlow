#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

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
  VFNEXT v1.3.1
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
  Run container:          ${ANSI_GREEN}${workflow.container}${ANSI_RESET}
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

  if (params.mode == "ILLUMINA"){
    ILLUMINA(reads_ch, ref_fa,ref_gff,ref_gcode)
    def bams_ch = ILLUMINA.out.align2ref_Out_ch.map{id, bams, _bais, _is_paired_end -> tuple(id, bams[0])}
    GENPLOTS(bams_ch)
  }

  if (params.mode == "NANOPORE"){
    NANOPORE(reads_ch, ref_fa)
    def bams_ch = NANOPORE.out.bam_ch.map{meta, bam -> tuple(meta.id, bam)}
    GENPLOTS(bams_ch)
  }

// -------------- Check if everything went okay -------------------------------
workflow.onComplete {
  if (workflow.success) {
    log.info """
      ===========================================
      ${ANSI_GREEN}Finished in ${workflow.duration}
      """.stripIndent()
      //See the report here ==> ${ANSI_RESET}$params.outDir/XXX_report.html

  } else {
    log.info """
      ===========================================
      ${ANSI_RED}Finished with errors!${ANSI_RESET}
      """.stripIndent()
  }
}

}


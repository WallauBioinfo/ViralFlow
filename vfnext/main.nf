#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import sub workflows
include { processInputs } from './workflows/step0-input-handling.nf'
include {ILLUMINA} from './workflows/ILLUMINA.nf'
include {NANOPORE} from './workflows/NANOPORE.nf'
include {GENPLOTS} from './workflows/GENPLOTS.nf'
include {METADATA} from './modules/metadata.nf'
include {
  metadataDir;
  metadataFailureMessage;
  writeRunManifest;
  containerSpecChannel;
  toolSpecChannel;
  referenceMetadataChannel
} from './modules/metadata_helpers.nf'

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

writeRunManifest(workflow, params, params.outDir, "RUNNING")

workflow.onError = {
  writeRunManifest(
    workflow,
    params,
    params.outDir,
    "FAILED",
    metadataFailureMessage(workflow)
  )
}

workflow.onComplete = {
  def finalStatus = workflow.success ? "SUCCESS" : "FAILED"
  writeRunManifest(
    workflow,
    params,
    params.outDir,
    finalStatus,
    workflow.success ? null : metadataFailureMessage(workflow)
  )

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

log.info """
  ===========================================
  VFNEXT ${workflow.manifest.version}
  parameters:
  -------------------------------------------
  --inDir            : ${params.inDir}
  --samplesheet      : ${params.samplesheet}
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
  readsCh = processInputs.out.readsCh
  refGff = processInputs.out.refGff
  refFa = processInputs.out.refFa
  refGcode = processInputs.out.refGcode

  readsMetadataCh = processInputs.out.sourceInputsCh

  referenceMetadataCh = referenceMetadataChannel("reference_fasta", refFa)

  gffMetadataCh = params.mode == "ILLUMINA"
    ? referenceMetadataChannel("referenceGff", refGff)
    : channel.empty()

  // Both modes clip primers when a BED is supplied, so both record it as a run
  // input.
  primerMetadataCh = params.primersBED
    ? channel.of(
        tuple(
          "reference",
          "primers_bed",
          file(params.primersBED).toAbsolutePath().normalize().toString(),
          file(params.primersBED)
        )
      )
    : channel.empty()

  samplesheetMetadataCh = params.samplesheet
    ? channel.of(
        tuple(
          "__run__",
          "samplesheet",
          file(params.samplesheet).toAbsolutePath().normalize().toString(),
          file(params.samplesheet)
        )
      )
    : channel.empty()

  checksumInputsCh = readsMetadataCh
    .concat(referenceMetadataCh)
    .concat(gffMetadataCh)
    .concat(primerMetadataCh)
    .concat(samplesheetMetadataCh)

  toolSpecsCh = toolSpecChannel(params, workflow, workflow.containerEngine)

  containerSpecsCh = containerSpecChannel(params, workflow.containerEngine)

  METADATA(
    checksumInputsCh,
    toolSpecsCh,
    containerSpecsCh,
    processInputs.out.resolvedInputsCh,
    metadataDir(params.outDir).toString()
  )

  if (params.mode == "ILLUMINA"){
    ILLUMINA(readsCh, refFa,refGff,refGcode)
    GENPLOTS(ILLUMINA.out.bamsCh, params.depth)
  }

  if (params.mode == "NANOPORE"){
    NANOPORE(readsCh, refFa)
    GENPLOTS(NANOPORE.out.bamsCh, params.np_min_depth)
  }

}

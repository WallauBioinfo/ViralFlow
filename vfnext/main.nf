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
  containerSpecs;
  toolSpecs
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

writeRunManifest(workflow, params, workflow.outputDir, "RUNNING")

workflow.onError = {
  writeRunManifest(
    workflow,
    params,
    workflow.outputDir,
    "FAILED",
    metadataFailureMessage(workflow)
  )
}

workflow.onComplete = {
  def finalStatus = workflow.success ? "SUCCESS" : "FAILED"
  writeRunManifest(
    workflow,
    params,
    workflow.outputDir,
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
  -output-dir        : ${workflow.outputDir}
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

  reads_metadata_ch = processInputs.out.source_inputs_ch

  reference_metadata_ch = ref_fa.map { reference ->
    tuple(
      "reference",
      "reference_fasta",
      reference.toAbsolutePath().normalize().toString(),
      reference
    )
  }

  gff_metadata_ch = params.mode == "ILLUMINA"
    ? ref_gff
        .filter { referenceGff -> referenceGff != null }
        .map { referenceGff ->
          tuple(
            "reference",
            "reference_gff",
            referenceGff.toAbsolutePath().normalize().toString(),
            referenceGff
          )
        }
    : channel.empty()

  primer_metadata_ch = params.primersBED
    ? channel.of(
        tuple(
          "reference",
          "primers_bed",
          file(params.primersBED).toAbsolutePath().normalize().toString(),
          file(params.primersBED)
        )
      )
    : channel.empty()

  samplesheet_metadata_ch = params.samplesheet
    ? channel.of(
        tuple(
          "__run__",
          "samplesheet",
          file(params.samplesheet).toAbsolutePath().normalize().toString(),
          file(params.samplesheet)
        )
      )
    : channel.empty()

  checksum_inputs_ch = reads_metadata_ch
    .concat(reference_metadata_ch)
    .concat(gff_metadata_ch)
    .concat(primer_metadata_ch)
    .concat(samplesheet_metadata_ch)

  tool_specs_ch = channel.fromList(
    toolSpecs(params, workflow).collect { spec ->
      tuple(spec.mode, spec.tool, spec.command, spec.container)
    }
  )

  container_specs_ch = channel.fromList(
    containerSpecs(params, workflow).collect { spec ->
      tuple(spec.name, spec.kind, spec.identity)
    }
  )

  METADATA(
    checksum_inputs_ch,
    tool_specs_ch,
    container_specs_ch,
    processInputs.out.resolved_inputs_ch,
    metadataDir(workflow.outputDir).toString()
  )

  if (params.mode == "ILLUMINA"){
    ILLUMINA(reads_ch, ref_fa,ref_gff,ref_gcode)
    GENPLOTS(ILLUMINA.out.bams_ch)
  }

  if (params.mode == "NANOPORE"){
    NANOPORE(reads_ch, ref_fa)
    GENPLOTS(NANOPORE.out.bams_ch)
  }

}

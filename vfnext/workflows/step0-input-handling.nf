#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {prepareDatabase} from "../modules/prepareDatabase.nf"

// set supported virus flag
def check_IL_custom_virus_params(errors) {
  def local_errors = errors

  // if a genome code was not provided, check if a gff and a ref fasta was
  if (params.refGenomeCode==null){
    if (params.runSnpEff==true){
      log.warn("The runSnpEff was set to ${params.runSnpEff}, but no refGenomeCode was provided.")
      log.warn("SnpEff will not be run")
    }
    if (params.referenceGFF==null){
      log.error("A 'custom' virus tag was set and no refGenomeCode was provided, therefore a referenceGFF must be provided.")
      local_errors += 1
    } else {
      def ref_gff_path = file(params.referenceGFF)
      if (!ref_gff_path.isFile()){
        log.error("${ref_gff_path} is not a file.")
        local_errors += 1
      }
      if (!ref_gff_path.exists()){
        log.error("${ref_gff_path} does not exist.")
        local_errors += 1
      }
    }

    if (params.referenceGenome==null){
      log.error("A 'custom' virus tag was set and no refGenomeCode was provided, therefore a referenceGenome must be provided.")
      local_errors += 1
    } else {
      def ref_fa_path = file(params.referenceGenome)
      if (!ref_fa_path.isFile()){
        log.error("${ref_fa_path} is not a file.")
        local_errors += 1
      }
      if (!ref_fa_path.exists()){
        log.error("${ref_fa_path} does not exists.")
        local_errors += 1
      }
    }
  }
  return local_errors
}

def validate_basic_params(accepted_modes) {
    def errors = 0
    
    if (!(params.mode in accepted_modes)) {
        log.error("The mode provided (${params.mode}) is not valid. Accepted modes are: ${accepted_modes.join(', ')}")
        errors += 1
    }
    
    return errors
}

def validate_primers_bed() {
  def errors = 0
  if (params.primersBED==null){
    //make adapter file optional, usefull for metagenomics
    log.warn("An BED file with primer positions was not provided. The pipeline will not run samtools clip to remove primer regions")
  }
  // if a path is provided, check if is valid
  else if (!(params.primersBED==null)){
    def adapter_fl = file(params.primersBED)
    if (!adapter_fl.isFile()){
      log.error("${params.primersBED} is not a file.")
      errors += 1
      }
    if (!adapter_fl.exists()){
      log.error("${params.primersBED} does not exists.")
      errors += 1
      }
  }
  return errors
}

def validate_virus_params() {
    def errors = 0
    def valid_virus = ["sars-cov2","custom"]
    
    if (!valid_virus.contains(params.virus)) {
        log.error("The virus provided (${params.virus}) is not valid.")
        errors += 1
    }

    // be sure custom only options were not set if a valid virus tag was provided
    if (valid_virus.contains(params.virus) && !(params.virus == "custom")) {
        if (!(params.referenceGFF==null)){
            log.warn("The valid virus tag (${params.virus}) was provided, ignoring the provided referenceGFF (${params.referenceGFF})")
            params.referenceGFF=null
        }
        if (!(params.referenceGenome==null)){
            log.warn("The valid virus tag (${params.virus}) was provided, ignoring the provided referenceGenome (${params.referenceGenome})")
            params.referenceGenome=null
        }
        if (!(params.refGenomeCode==null)){
            log.warn("The valid virus tag (${params.virus}) was provided, ignoring the provided refGenomeCode (${params.refGenomeCode})")
            params.refGenomeCode=null
        }
    }

    // if a custom virus, check if mandatory params were set
    if (params.virus=="custom"){
        errors += check_IL_custom_virus_params(errors)
    }
    
    return errors
}

def validate_illumina_params() {
    def errors = 0
    
    // Primer BED validation
    errors += validate_primers_bed()
    
    // Virus validation
    errors += validate_virus_params()
    
    return errors
}

def validate_nanopore_params() {
    def errors = 0
    
    if (!params.referenceGenome) {
        log.error("A reference genome fasta file must be provided for NANOPORE mode")
        errors += 1
    } else {
        def ref_fa_path = file(params.referenceGenome)
        if (!ref_fa_path.exists() || !ref_fa_path.isFile()) {
            log.error("Reference genome file ${params.referenceGenome} does not exist or is not a file")
            errors += 1
        }
    }
    
    return errors
}
def validate_directories() {
  def errors = 0
  
  // check if output dir exists, if not create the default
  if (params.outDir){
    def outDir_path = file(params.outDir)
  
    if (!outDir_path.exists()){
      log.warn("${params.outDir} does not exist, the directory will be created")
      outDir_path.mkdir()
    }
    if (!(outDir_path.isDirectory())){
      log.error("${params.outDir} is not a directory")
      errors+=1
    }
  }

  // check if input dir exists
  if (params.inDir==null){
    log.error("An input directory must be provided.")
    errors+=1
  }
  if (params.inDir){
    def inDir_path = file(params.inDir)
    if (!inDir_path.isDirectory()){
      log.error("${params.inDir} is not a directory")
      errors+=1
    }
  }
  return errors

}

def validate_resources() {
  def errors = 0
  // get number of cpus available for nextflow if running local
  def maxcpus = Runtime.runtime.availableProcessors()

  if (workflow.profile == "standard"){
    // if cpus were not specified or higher than the available cpus, set it to use all cpus available
    if ((params.nextflowSimCalls == null) || (params.nextflowSimCalls > maxcpus)){
      log.warn("Number of requested simultaneous nextflow calls (${params.nextflowSimCalls}) was set to max cpus available (${maxcpus})")
      params.nextflowSimCalls = maxcpus
    }
  }

  // top multithread to maxcpus
  if (params.fastp_threads > maxcpus){
    log.warn("Number of threads to be used by fastp (${params.fastp_threads}) is higher than requested cpus (${maxcpus}). Setting it to ${maxcpus}.")
    params.fastp_threads = maxcpus
  }

  if (params.bwa_threads > maxcpus){
    log.warn("Number of threads to be used by bwa (${params.bwa_threads}) is higher than available threads (${maxcpus}). Setting it to ${maxcpus}.")
    params.bwa_threads = maxcpus
  }

  if (params.mafft_threads > maxcpus){
    log.warn("Number of threads to be used by mafft (${params.mafft_threads}) is higher than available threads (${maxcpus}). Setting it to ${maxcpus}.")
    params.mafft_threads = maxcpus
  }
    return errors
}

def validate_parameters() {
    def errors = 0
    def ACCEPTED_MODES = ["ILLUMINA", "NANOPORE"]
    
    // Basic parameter validation
    errors += validate_basic_params(ACCEPTED_MODES)
    
    // Mode-specific validation
    if (params.mode == "ILLUMINA") {
        errors += validate_illumina_params()
    } else if (params.mode == "NANOPORE") {
        errors += validate_nanopore_params()
    }
    
    // Common validation
    errors += validate_directories()
    errors += validate_resources()
    
    // Exit if errors found
    if (errors > 0) {
        log.error("${errors} validation errors detected")
        exit 1
    }
}

workflow processInputs {
  main:
    // --- Sanity Check -------------------------------------------------------
    // check if fasta exists and follow symlinks if needed
    //ref_fa = file(reference_fasta, checkIfExists=true, followLinks=true)
    //-------------------------------------------------------------------------
    validate_parameters()
    if (params.mode == "ILLUMINA"){
      // ---- get reference GFF and fasta ---------------------------------------
      // Setup ref code values for supported virus
      ref_gcode = null
      reference_fa = null
      reference_gff = null

      if (!(params.virus=="custom")){
        if (params.virus=="sars-cov2"){
          ref_gcode = "NC_045512.2"
        }
      }

      // if custom virus, check if a genome code was provided, if not
      // emit the ref gff and fasta provided
      if (params.virus=="custom"){
        if (!(params.refGenomeCode==null)){
          ref_gcode = params.refGenomeCode
        } else {
          reference_gff = params.referenceGFF
          reference_fa = params.referenceGenome
        }
      }

      // if a genome code was provided, get the reference fasta and gff
      if (!(ref_gcode==null)){
        prepareDatabase(ref_gcode)
        reference_fa = prepareDatabase.out.ref_fa
        reference_gff = prepareDatabase.out.ref_gff
      }

      // be sure a reference fasta and a reference gff was obtained
      assert !(reference_fa == null) && !(reference_gff == null)
    }

    if (params.mode == "NANOPORE"){
      // if a reference fasta was provided, use it
      if (params.referenceGenome){
        reference_fa = file(params.referenceGenome)
      } else {
        log.error("A reference genome fasta file must be provided for NANOPORE mode")
        exit 1
      }

      reference_gff = null
      ref_gcode = null
    }
    // get reads
    // current support follow the rules:
    // paired reads with R1 AND R2 pattern and .fq.gz / .fastq.gz extensions
    // single reads for files without R1 or R2 pattern and .fq.gz / .fastq.gz extensions
    // if some file has 0 bytes, the file is removed of the analysis.
    
    reads_channel_paired_raw = Channel
      .fromFilePairs(["${params.inDir}/*_R{1,2}*.fq.gz", "${params.inDir}/*_R{1,2}*.fastq.gz"])  
    reads_channel_paired_raw
      .filter{(it[1][0].size()==0) && (it[1][1].size()==0)}
      .view{log.warn("Excluding ${it[0]} fastq files due to 0 bytes size")}    
    reads_channel_paired = reads_channel_paired_raw.filter{(it[1][0].size()>0) && (it[1][1].size()>0)}

    reads_channel_single_raw = channel
      .fromPath(["${params.inDir}/*.fq.gz", "${params.inDir}/*.fastq.gz", "${params.inDir}/*.fastq"])
      .collect()
      .map { files ->
          def grouped = files.groupBy { file ->
              file.getName().replaceAll('_R[12]', '')
          }

          grouped.collectMany { _sample, fileList ->
              def hasR2 = fileList.any { it.getName().contains('_R2') }
              hasR2 ? fileList.findAll { !it.getName().contains('_R1') && !it.getName().contains('_R2') } : fileList
          }
      }
      .flatten()
      .map { file -> 
          def fileName = file.getName()
          def baseName = fileName.contains('_R1') ? fileName.split('_R1')[0] : fileName.split('\\.')[0]
          [baseName, [file]] 
      }
    reads_channel_single_raw
      .filter{it[1][0].size() == 0}
      .view{log.warn("Excluding ${it[0]} fastq files due to 0 bytes size")}
    reads_channel_single = reads_channel_single_raw.filter{(it[1][0].size()>0)}

    reads_channel = reads_channel_paired
                    .concat(reads_channel_single)
                    .map { sample_id, files -> 
                          def meta = [id: sample_id,
                                      is_paired_end: files.size() == 2]
                          tuple(meta, files)
                    }
  emit:
    reads_ch = reads_channel
    ref_gff = reference_gff
    ref_fa = reference_fa
    ref_gcode = ref_gcode
}

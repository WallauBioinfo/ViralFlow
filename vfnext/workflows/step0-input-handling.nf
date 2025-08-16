#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {prepareDatabase} from "../modules/prepareDatabase.nf"

// set supported virus flag

def validate_parameters() {
    // --- SANITY CHECKS ------------------------------------------------------
    def errors = 0
    // check if required params were provided and if files provided exists

    if (params.primersBED==null){
      //make adapter file optional, usefull for metagenomics
      log.warn("An BED file with primer positions was not provided. The pipeline will not run samtools clip to remove primer regions")
      }
    // if only the flag is provided withou any value, it is considered as true
    else if (params.primersBED==true){
      log.error("the BED file flag was set but no value provided")
      errors +=1
    }
    // if a path is provided, check if is valid
    else if (!(params.primersBED==null)){
        def adapter_fl = file(params.primersBED)
        if (!adapter_fl.isFile()){
          log.error("${params.primersBED} is not a file.")
          errors += 1
        }
      //errors +=1
    }
    // --- VIRUS FLAGS CHECK
    // check if a valid virus flag was provided
    def valid_virus = ["sars-cov2","custom"]
    if (!valid_virus.contains(params.virus)) {
      log.error("The virus provided (${params.virus}) is not valid.")
      errors += 1
    }

    // be sure custom only options were not set if a valid virus tag was provided
    if (valid_virus.contains(params.virus) && !(params.virus == "custom")) {
        if (!(params.referenceGFF==null)){
          log.warn("The valid virus tag (${params.virus}) was provided, ingnoring the provided referenceGFF (${params.referenceGFF})")
          params.referenceGFF=null
        }
        if (!(params.referenceGenome==null)){
          log.warn("The valid virus tag (${params.virus}) was provided, ingnoring the provided referenceGenome (${params.referenceGenome})")
          params.referenceGenome=null
        }
        if (!(params.refGenomeCode==null)){
          log.warn("The valid virus tag (${params.virus}) was provided, ingnoring the provided refGenomeCode (${params.refGenomeCode})")
          params.refGenomeCode=null
        }

    }
    // ------------------------------------------------------------------------
    // if a custom virus, check if mandatory params were set
    if (params.virus=="custom"){
      // if a genome code was not provided, check if a gff and a ref fasta was
      if (params.refGenomeCode==null){
        if (params.runSnpEff==true){
          log.warn("The runSnpEff was set to ${params.runSnpEff}, but no refGenomeCode was provided.")
          log.warn("SnpEff will not be run")
        }
        if (params.referenceGFF==null){
          log.error("A 'custom' virus tag was set and no refGenomeCode was provided, therefore a referenceGFF must be provided.")
          errors += 1
        } else {
          def ref_gff_path = file(params.referenceGFF)
          if (!ref_gff_path.isFile()){
            log.error("${ref_gff_path} is not a file.")
            errors += 1
          }
          if (!ref_gff_path.exists()){
            log.error("${ref_gff_path} does not exists.")
            errors += 1
          }
        }

        if (params.referenceGenome==null){
          log.error("A 'custom' virus tag was set and no refGenomeCode was provided, therefore a referenceGenome must be provided.")
          errors += 1
        } else {
          def ref_fa_path = file(params.referenceGenome)
          if (!ref_fa_path.isFile()){
            log.error("${ref_fa_path} is not a file.")
            errors += 1
          }
          if (!ref_fa_path.exists()){
            log.error("${ref_fa_path} does not exists.")
            errors += 1
          }
        }
      }
    }
    // ------------------------------------------------------------------------

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

    //TODO check if adapterFile exist
    //-------------------------------------------------------------------------
    // count errors and kill nextflow if any had been found
    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
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

    // get reads
    // current support follow the rules:
    // paired reads with R1 AND R2 pattern and .fq.gz / .fastq.gz extensions
    // single reads for files without R1 or R2 pattern and .fq.gz / .fastq.gz extensions
    // if some file has 0 bytes, the file is removed of the analysis.

    reads_channel_paired_raw = channel
      .fromFilePairs(["${params.inDir}/*_R{1,2}*.fq.gz", "${params.inDir}/*_R{1,2}*.fastq.gz"])  
    
    reads_channel_paired_raw
      .filter((it[1][0].size()==0) && (it[1][1].size()==0))
      .view(log.warn("Excluding ${it[0]} fastq files due to 0 bytes size"))    
    reads_channel_paired = reads_channel_paired_raw.filter((it[1][0].size()>0) && (it[1][1].size()>0))

    reads_channel_single_raw = channel
      .fromPath(["${params.inDir}/*.fq.gz", "${params.inDir}/*.fastq.gz"])
      .collect()
      .map { files ->
          def grouped = files.groupBy { file ->
              file.getName().replaceAll('_R[12]', '')
          }

          grouped.collectMany { sample, fileList ->
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
      .filter(it[1][0].size() == 0)
      .view(log.warn("Excluding ${it[0]} fastq files due to 0 bytes size"))
    reads_channel_single = reads_channel_single_raw.filter((it[1][0].size()>0))

    reads_channel = reads_channel_paired.concat(reads_channel_single)

  emit:
    reads_ch = reads_channel
    ref_gff = reference_gff
    ref_fa = reference_fa
    ref_gcode = ref_gcode
}

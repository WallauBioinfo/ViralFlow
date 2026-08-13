#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {prepareDatabase} from "../modules/prepareDatabase.nf"

def parseCsvRecord(String line) {
  if ((line.count('"') % 2) != 0) throw new IllegalArgumentException('unterminated quoted field')
  line.split(/,(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)/, -1).toList().collect { raw ->
    def value = raw.trim()
    if (value.startsWith('"') && value.endsWith('"') && value.size() >= 2) {
      value = value.substring(1, value.size() - 1).replace('""', '"')
    } else if (value.contains('"')) {
      throw new IllegalArgumentException("invalid quote placement in '${value}'")
    }
    value
  }
}

def canonicalPath(java.nio.file.Path baseDir, String rawPath) {
  def candidate = java.nio.file.Path.of(rawPath)
  if (!candidate.isAbsolute()) candidate = baseDir.resolve(candidate)
  candidate.toAbsolutePath().normalize()
}

def validateFastqPath(java.nio.file.Path path, String location, List errors) {
  if (!java.nio.file.Files.exists(path)) {
    errors << "${location}: FASTQ does not exist: ${path}"
    return null
  }
  if (!java.nio.file.Files.isRegularFile(path)) {
    errors << "${location}: FASTQ is not a regular file: ${path}"
    return null
  }
  if (!java.nio.file.Files.isReadable(path)) errors << "${location}: FASTQ is not readable: ${path}"
  if (java.nio.file.Files.size(path) == 0) errors << "${location}: FASTQ is empty: ${path}"
  if (!(path.fileName.toString() ==~ /(?i).+\.(fastq|fq)(\.gz)?/)) {
    errors << "${location}: unsupported FASTQ extension: ${path}"
  }
  try {
    return path.toRealPath()
  } catch (Exception _ignored) {
    return path
  }
}

def validateCanonicalRows(List rows, String mode) {
  def errors = []
  def safeId = ~/[A-Za-z0-9][A-Za-z0-9._-]*/
  def safeColumn = ~/[A-Za-z_][A-Za-z0-9_]*/
  def reserved = ['id', 'is_paired_end', 'chunk_index'] as Set
  def usedPaths = [:]
  def sampleLayouts = [:].withDefault { [] }
  def sampleMetadata = [:]

  rows.each { row ->
    def location = row.location
    if (!row.sample_id || !(row.sample_id ==~ safeId)) {
      errors << "${location}: unsafe sample_id '${row.sample_id ?: ''}'; expected [A-Za-z0-9][A-Za-z0-9._-]*"
    }
    row.metadata.each { key, _value ->
      if (!(key ==~ safeColumn) || reserved.contains(key)) {
        errors << "${location}: invalid or reserved metadata column '${key}'"
      }
    }
    def path1 = validateFastqPath(row.fastq_1, "${location} fastq_1", errors)
    def path2 = row.fastq_2 ? validateFastqPath(row.fastq_2, "${location} fastq_2", errors) : null
    row.fastq_1 = path1 ?: row.fastq_1
    row.fastq_2 = path2

    [fastq_1: path1, fastq_2: path2].each { role, path ->
      if (path) {
        def identity = path.toString()
        if (usedPaths.containsKey(identity)) {
          errors << "${location}: ${path} is already assigned at ${usedPaths[identity]}"
        } else {
          usedPaths[identity] = "${location} ${role}"
        }
      }
    }

    def paired = path2 != null
    sampleLayouts[row.sample_id] << paired
    if (mode == 'NANOPORE' && paired) errors << "${location}: fastq_2 must be empty in NANOPORE mode"

    if (!sampleMetadata.containsKey(row.sample_id)) {
      sampleMetadata[row.sample_id] = row.metadata
    } else if (sampleMetadata[row.sample_id] != row.metadata) {
      errors << "${location}: metadata values are inconsistent for sample '${row.sample_id}'"
    }
  }

  sampleLayouts.each { sampleId, layouts ->
    if (layouts.unique().size() > 1) errors << "sample '${sampleId}' mixes single-end and paired-end rows"
  }
  if (!rows) errors << 'No FASTQ inputs were found'
  if (errors) error "Input validation failed:\n - ${errors.join('\n - ')}"
  return rows
}

def parseSamplesheet(String samplesheet, String mode) {
  def sheetPath = java.nio.file.Path.of(samplesheet).toAbsolutePath().normalize()
  if (!java.nio.file.Files.exists(sheetPath) || !java.nio.file.Files.isRegularFile(sheetPath)) {
    error "Sample sheet does not exist or is not a file: ${sheetPath}"
  }
  def lines = java.nio.file.Files.readAllLines(sheetPath)
  if (!lines) error "Sample sheet is empty: ${sheetPath}"
  if (lines[0].startsWith('\uFEFF')) lines[0] = lines[0].substring(1)
  def headers
  try { headers = parseCsvRecord(lines[0]) }
  catch (Exception exception) { error "Invalid sample-sheet header: ${exception.message}" }
  if (headers.size() != headers.unique().size()) error 'Sample sheet contains duplicate column names'
  def required = ['sample_id', 'fastq_1', 'fastq_2']
  def missing = required.findAll { header -> !headers.contains(header) }
  if (missing) error "Sample sheet is missing required columns: ${missing.join(', ')}"
  def metadataHeaders = headers.findAll { header -> !required.contains(header) }
  def baseDir = sheetPath.parent
  def rows = []
  lines.drop(1).eachWithIndex { line, index ->
    if (!line.trim()) return
    def values
    try { values = parseCsvRecord(line) }
    catch (Exception exception) { error "Sample sheet row ${index + 2}: ${exception.message}" }
    if (values.size() != headers.size()) error "Sample sheet row ${index + 2} has ${values.size()} fields; expected ${headers.size()}"
    def record = [headers, values].transpose().collectEntries()
    if (!record.fastq_1) error "Sample sheet row ${index + 2}: fastq_1 is required"
    rows << [
      sample_id: record.sample_id,
      fastq_1: canonicalPath(baseDir, record.fastq_1),
      fastq_2: record.fastq_2 ? canonicalPath(baseDir, record.fastq_2) : null,
      metadata: metadataHeaders.collectEntries { header -> [(header): record[header]] },
      location: "sample-sheet row ${index + 2}"
    ]
  }
  validateCanonicalRows(rows, mode)
}

def parseLegacyDirectory(String inDir, String mode) {
  log.warn("--inDir automatic discovery is deprecated and will be removed in ViralFlow v3; use --samplesheet")
  def inputPath = java.nio.file.Path.of(inDir).toAbsolutePath().normalize()
  if (!java.nio.file.Files.isDirectory(inputPath)) error "${inputPath} is not a directory"
  def files = []
  java.nio.file.Files.list(inputPath).withCloseable { stream ->
    stream.filter { path -> java.nio.file.Files.isRegularFile(path) }
      .filter { path -> path.fileName.toString() ==~ /(?i).+\.(fastq|fq)(\.gz)?/ }
      .sorted()
      .forEach { path -> files << path }
  }
  def paired = [:].withDefault { [:] }
  def singles = []
  files.each { path ->
    def matcher = path.fileName.toString() =~ /(?i)^(.+)_R([12])(?:_[^.]+)?\.(fastq|fq)(\.gz)?$/
    if (matcher.matches()) {
      def sampleId = matcher[0][1]
      def mate = matcher[0][2]
      if (paired[sampleId].containsKey(mate)) {
        error "Legacy input discovery produced duplicate sample ID '${sampleId}' for mate R${mate}; use --samplesheet to define chunks explicitly"
      }
      paired[sampleId][mate] = path
    }
    else singles << path
  }
  def rows = []
  paired.each { sampleId, mates ->
    if (!(mates['1'] && mates['2'])) error "Legacy input sample '${sampleId}' has an orphan Illumina mate"
    rows << [sample_id: sampleId, fastq_1: mates['1'], fastq_2: mates['2'], metadata: [:], location: "legacy sample ${sampleId}"]
  }
  singles.each { path ->
    def sampleId = path.fileName.toString().replaceFirst(/(?i)\.(fastq|fq)(\.gz)?$/, '')
    rows << [sample_id: sampleId, fastq_1: path, fastq_2: null, metadata: [:], location: "legacy file ${path.fileName}"]
  }
  def duplicateIds = rows.groupBy { row -> row.sample_id }.findAll { _sample_id, sample_rows -> sample_rows.size() > 1 }.keySet()
  if (duplicateIds) error "Legacy input discovery produced duplicate sample IDs: ${duplicateIds.sort().join(', ')}; use --samplesheet to define chunks explicitly"
  validateCanonicalRows(rows.sort { row -> row.sample_id }, mode)
}

def groupCanonicalRows(List rows) {
  def grouped = new LinkedHashMap()
  rows.each { row ->
    row.chunk_index = (grouped[row.sample_id]?.size() ?: 0) + 1
    grouped.computeIfAbsent(row.sample_id) { [] } << row
  }
  grouped.collect { sampleId, chunks ->
    def paired = chunks[0].fastq_2 != null
    def meta = [id: sampleId, is_paired_end: paired] + chunks[0].metadata
    tuple(
      meta,
      chunks.collect { chunk -> file(chunk.fastq_1.toString()) },
      paired ? chunks.collect { chunk -> file(chunk.fastq_2.toString()) } : []
    )
  }
}

process prepare_sample_reads {
  tag "${meta.id}"

  input:
    tuple val(meta), path(fastq_1_chunks, stageAs: 'r1/chunk??/*'), path(fastq_2_chunks, stageAs: 'r2/chunk??/*')

  output:
    tuple val(meta), path("${meta.id}.*.fastq.gz")

  script:
    def r1Inputs = fastq_1_chunks instanceof List ? fastq_1_chunks : [fastq_1_chunks]
    def r2Inputs = fastq_2_chunks instanceof List ? fastq_2_chunks : (fastq_2_chunks ? [fastq_2_chunks] : [])
    def r1 = r1Inputs.collect { input ->
      input.name.toLowerCase().endsWith('.gz') ? "gzip -cd '${input}'" : "cat '${input}'"
    }.join('; ')
    def r2 = r2Inputs.collect { input ->
      input.name.toLowerCase().endsWith('.gz') ? "gzip -cd '${input}'" : "cat '${input}'"
    }.join('; ')
    if (meta.is_paired_end) {
      """
      set -euo pipefail
      { ${r1}; } | gzip -n -c > '${meta.id}.R1.fastq.gz'
      { ${r2}; } | gzip -n -c > '${meta.id}.R2.fastq.gz'
      gzip -t '${meta.id}.R1.fastq.gz' '${meta.id}.R2.fastq.gz'
      """
    } else {
      """
      set -euo pipefail
      { ${r1}; } | gzip -n -c > '${meta.id}.SE.fastq.gz'
      gzip -t '${meta.id}.SE.fastq.gz'
      """
    }
}

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

  def configuredOutputDir = java.nio.file.Path.of(params.outDir.toString()).toAbsolutePath().normalize()
  def runtimeOutputDir = java.nio.file.Path.of(workflow.outputDir.toString()).toAbsolutePath().normalize()
  if (configuredOutputDir != runtimeOutputDir) {
    log.error("Conflicting output directories detected: use --outDir instead of Nextflow -output-dir or an outputDir config override")
    errors += 1
  }

  // check if output dir exists, if not create the default
  if (params.outDir){
    def outDir_path = file(params.outDir)

    if (!outDir_path.exists()){
      log.warn("${params.outDir} does not exist, the directory will be created")
      outDir_path.mkdirs()
    }
    if (!(outDir_path.isDirectory())){
      log.error("${params.outDir} is not a directory")
      errors+=1
    }
  }

  if (params.samplesheet && params.inDir) {
    log.error("--samplesheet and --inDir cannot be used together")
    errors += 1
  }
  if (params.samplesheet) {
    def sheetPath = file(params.samplesheet)
    if (!sheetPath.isFile()) {
      log.error("${params.samplesheet} is not a sample-sheet file")
      errors += 1
    }
  } else if (params.inDir) {
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
        error "${errors} validation errors detected"
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
        error "A reference genome fasta file must be provided for NANOPORE mode"
      }

      reference_gff = null
      ref_gcode = null
    }
    def effectiveInDir = params.inDir ?: workflow.launchDir.resolve('input').toString()
    def canonicalRows = params.samplesheet
      ? parseSamplesheet(params.samplesheet.toString(), params.mode.toString())
      : parseLegacyDirectory(effectiveInDir, params.mode.toString())
    def groupedInputs = groupCanonicalRows(canonicalRows)
    prepare_sample_reads(channel.fromList(groupedInputs))
    prepared_reads = prepare_sample_reads.out.map { meta, reads ->
      tuple(meta, reads instanceof List ? reads : [reads])
    }

    source_inputs = channel.fromList(canonicalRows.collectMany { row ->
      def records = [tuple(row.sample_id, "fastq_1_chunk_${row.chunk_index}", row.fastq_1.toString(), file(row.fastq_1.toString()))]
      if (row.fastq_2) records << tuple(row.sample_id, "fastq_2_chunk_${row.chunk_index}", row.fastq_2.toString(), file(row.fastq_2.toString()))
      records
    })

    resolved_inputs = channel.of(canonicalRows.collect { row ->
      [row.sample_id, row.chunk_index, row.fastq_2 ? 'paired' : 'single', row.fastq_1.toString(), row.fastq_2?.toString() ?: '']
    })
  emit:
    reads_ch = prepared_reads
    source_inputs_ch = source_inputs
    resolved_inputs_ch = resolved_inputs
    ref_gff = reference_gff
    ref_fa = reference_fa
    ref_gcode = ref_gcode
}

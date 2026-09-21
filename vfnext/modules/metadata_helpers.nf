def metadataDir(outputDir) {
    java.nio.file.Path.of(outputDir.toString()).toAbsolutePath().normalize()
        .resolve('RUN_METADATA')
}

def normalizeMetadata(value) {
    if (value == null || value instanceof Boolean || value instanceof Number || value instanceof String) {
        return value
    }
    if (value instanceof Path || value instanceof File) {
        return absoluteMetadataPath(value)
    }
    if (value instanceof Date) {
        return value.toInstant().toString()
    }
    if (value instanceof Map) {
        return value.collectEntries { key, item -> [(key.toString()): normalizeMetadata(item)] }
    }
    if (value instanceof Collection) {
        return value.collect { item -> normalizeMetadata(item) }
    }
    if (value.getClass().isArray()) {
        return value.toList().collect { item -> normalizeMetadata(item) }
    }
    return value.toString()
}

def absoluteMetadataPath(value) {
    if (value == null) {
        return null
    }
    try {
        return java.nio.file.Path.of(value.toString()).toAbsolutePath().normalize().toString()
    }
    catch (Exception _ignored) {
        return value.toString()
    }
}

def safeMetadataValue(closure) {
    try {
        return closure.call()
    }
    catch (Exception _ignored) {
        return null
    }
}

def gitMetadataValue(projectDir, arguments) {
    try {
        def command = ['git', '-C', projectDir.toString()] + arguments
        def process = new ProcessBuilder(command).redirectErrorStream(true).start()
        def output = process.inputStream.text.trim()
        return process.waitFor() == 0 ? output : null
    }
    catch (Exception _ignored) {
        return null
    }
}

def metadataFailureMessage(workflow) {
    safeMetadataValue { -> workflow.errorMessage } ?: safeMetadataValue { -> workflow.errorReport }
}

def writeRunManifest(workflow, params, configuredOutputDir, status, failureMessage = null) {
    def metadataOutputDir = metadataDir(configuredOutputDir)
    java.nio.file.Files.createDirectories(metadataOutputDir)

    def profile = safeMetadataValue { -> workflow.profile }?.toString() ?: ''
    def gitStatus = gitMetadataValue(workflow.projectDir, ['status', '--porcelain'])
    def manifest = [
        schema_version: 1,
        pipeline: [
            name: safeMetadataValue { -> workflow.manifest.name },
            version: safeMetadataValue { -> workflow.manifest.version },
            repository: safeMetadataValue { -> workflow.manifest.homePage },
            revision: gitMetadataValue(workflow.projectDir, ['rev-parse', '--abbrev-ref', 'HEAD']),
            commit_id: gitMetadataValue(workflow.projectDir, ['rev-parse', 'HEAD']),
            git_dirty: gitStatus == null ? null : !gitStatus.isEmpty()
        ],
        execution: [
            status: status,
            failure_message: failureMessage,
            session_id: safeMetadataValue { -> workflow.sessionId },
            run_name: safeMetadataValue { -> workflow.runName },
            command_line: safeMetadataValue { -> workflow.commandLine },
            profile: profile,
            nextflow_version: nextflow.BuildInfo.version,
            start_time: normalizeMetadata(safeMetadataValue { -> workflow.start }),
            end_time: status == 'RUNNING' ? null : normalizeMetadata(new Date()),
            duration: normalizeMetadata(safeMetadataValue { -> workflow.duration }),
            success: status == 'SUCCESS'
        ],
        runtime: [
            user: safeMetadataValue { -> workflow.userName },
            host: safeMetadataValue { -> java.net.InetAddress.localHost.hostName },
            os: System.getProperty('os.name'),
            os_version: System.getProperty('os.version'),
            architecture: System.getProperty('os.arch'),
            executor: profile.contains('pbs') ? 'pbs' : 'local',
            container_engine: profile.contains('apptainer') ? 'apptainer' : 'singularity'
        ],
        paths: [
            launch_dir: absoluteMetadataPath(safeMetadataValue { -> workflow.launchDir }),
            project_dir: absoluteMetadataPath(safeMetadataValue { -> workflow.projectDir }),
            work_dir: absoluteMetadataPath(safeMetadataValue { -> workflow.workDir }),
            input_dir: params.samplesheet
                ? null
                : absoluteMetadataPath(params.inDir ?: workflow.launchDir.resolve('input')),
            samplesheet: absoluteMetadataPath(params.samplesheet),
            output_dir: absoluteMetadataPath(configuredOutputDir)
        ],
        analysis: [
            mode: normalizeMetadata(params.mode),
            virus: normalizeMetadata(params.virus),
            clair3_model: params.mode == 'NANOPORE' ? normalizeMetadata(params.clair3_model) : null,
            clair3_qual: params.mode == 'NANOPORE' ? normalizeMetadata(params.clair3_qual) : null,
            mapping_quality: params.mode == 'NANOPORE' ? normalizeMetadata(params.mapping_quality) : null,
            af_threshold: params.mode == 'NANOPORE' ? normalizeMetadata(params.af_threshold) : null,
            min_depth: params.mode == 'NANOPORE' ? normalizeMetadata(params.np_min_depth) : normalizeMetadata(params.depth),
            consensus_mask_rule: params.mode == 'NANOPORE' ? 'depth <= min_depth' : null
        ],
        parameters: normalizeMetadata(params.entrySet().collectEntries { entry ->
            [(entry.key.toString()): entry.value]
        }),
        files: [
            input_checksums: 'input_checksums.tsv',
            resolved_sample_inputs: 'resolved_sample_inputs.tsv',
            software_versions: 'software_versions.tsv',
            containers: 'container_manifest.tsv',
            trace: metadataOutputDir.resolve('execution_trace.tsv').toString(),
            report: metadataOutputDir.resolve('execution_report.html').toString(),
            timeline: metadataOutputDir.resolve('execution_timeline.html').toString()
        ]
    ]

    def target = metadataOutputDir.resolve('run_manifest.json')
    def temporary = metadataOutputDir.resolve('run_manifest.json.tmp')
    temporary.toFile().text = groovy.json.JsonOutput.prettyPrint(
        groovy.json.JsonOutput.toJson(manifest)
    ) + System.lineSeparator()
    java.nio.file.Files.move(
        temporary,
        target,
        java.nio.file.StandardCopyOption.REPLACE_EXISTING,
        java.nio.file.StandardCopyOption.ATOMIC_MOVE
    )
}

// Build checksum records for a reference input channel. Shared by main.nf and
// the metadata tests so the tests exercise the production path rather than a
// copy of it. file() normalises the value: callers may hand over a Path or a
// raw parameter string, and null entries are dropped for optional inputs.
def referenceMetadataChannel(role, source) {
    source
        .filter { value -> value != null }
        .map { value ->
            def resolved = file(value)
            tuple(
                "reference",
                role,
                resolved.toAbsolutePath().normalize().toString(),
                resolved
            )
        }
}

def localContainerSpec(name, pathValue) {
    def identity = absoluteMetadataPath(pathValue)
    def kind = java.nio.file.Files.isDirectory(java.nio.file.Path.of(identity)) ? 'local_sandbox' : 'local_sif'
    [name: name, kind: kind, identity: identity]
}

def remoteContainerSpec(name, uri) {
    [name: name, kind: 'remote_uri', identity: uri]
}

def toolSpec(mode, name, command, containerValue) {
    [mode: mode, tool: name, command: command, container: containerValue.toString()]
}

// The channel builders below own the spec-map -> tuple mapping. capture_tool_version
// and capture_container_metadata read those tuples positionally, so keeping the
// mapping in one place is what lets a test pin the field order.
def containerSpecChannel(params, workflow) {
    channel.fromList(
        containerSpecs(params, workflow).collect { spec ->
            tuple(spec.name, spec.kind, spec.identity)
        }
    )
}

def toolSpecChannel(params, workflow) {
    channel.fromList(
        toolSpecs(params, workflow).collect { spec ->
            tuple(spec.mode, spec.tool, spec.command, spec.container)
        }
    )
}

// Resolves an image by the name configs/containers.config declares it under.
// Failing loudly on an unknown name is the point: the alternative is a manifest
// that omits an image, or names one the run never used.
def illuminaContainerSpec(params, name) {
    def image = (params.illumina_containers ?: [:])[name]
    if (!image) {
        error "No image declared for '${name}' in params.illumina_containers (configs/containers.config)"
    }
    return localContainerSpec(name, image)
}

// Takes workflow to match toolSpecs() and containerSpecChannel(), though it no
// longer needs it: the ILLUMINA paths used to be built from workflow.projectDir
// and now come from params.illumina_containers.
def containerSpecs(params, _workflow) {
    def specs = []
    if (params.mode == 'NANOPORE') {
        specs << localContainerSpec('nanopore_base', params.base_container)
        specs << remoteContainerSpec('clair3', params.clair3_container)
    }
    else if (params.mode == 'ILLUMINA') {
        // Names, not paths: the paths live in params.illumina_containers, which
        // the process directives read too, so the manifest cannot drift from
        // the images the run actually used. What stays here is which of those
        // images a given run reaches for, since that depends on the parameters
        // below rather than on the configuration.
        def names = [
            'edirect',
            'generate_consensus',
            'fastp',
            'samtools',
            'mafft',
            'picard',
            'intrahost_analysis',
            'generate_plots',
            'compiled_outputs'
        ]
        if (params.runSnpEff) {
            names += ['snpeff', 'generate_report']
        }
        if (params.virus == 'sars-cov2') {
            names += ['pangolin', 'nextclade']
        }
        specs.addAll(names.collect { name -> illuminaContainerSpec(params, name) })
    }
    specs.unique { spec -> spec.identity }
}

def toolSpecs(params, workflow) {
    if (params.mode == 'NANOPORE') {
        return [
            toolSpec('NANOPORE', 'porechop_abi', 'porechop_abi --version', params.base_container),
            toolSpec('NANOPORE', 'minimap2', 'minimap2 --version', params.base_container),
            toolSpec('NANOPORE', 'samtools', 'samtools --version | head -n 1', params.base_container),
            toolSpec('NANOPORE', 'bcftools', 'bcftools --version | head -n 1', params.base_container),
            toolSpec('NANOPORE', 'clair3', 'run_clair3.sh -v ', params.clair3_container)
        ]
    }

    def containerDir = java.nio.file.Path.of(workflow.projectDir.toString()).resolve('containers')
    def specs = [
        toolSpec('ILLUMINA', 'fastp', 'fastp --version', containerDir.resolve('fastp:1.0.1.sif')),
        toolSpec('ILLUMINA', 'bwa', 'bwa 2>&1 | head -n 3', containerDir.resolve('generate_consensus:2.0.0.sif')),
        toolSpec('ILLUMINA', 'samtools', 'samtools --version | head -n 1', containerDir.resolve('generate_consensus:2.0.0.sif')),
        toolSpec('ILLUMINA', 'ivar', 'ivar version', containerDir.resolve('generate_consensus:2.0.0.sif')),
        toolSpec('ILLUMINA', 'mafft', 'mafft --version', containerDir.resolve('mafft:7.505_2.sif'))
    ]
    if (params.runSnpEff) {
        specs << toolSpec('ILLUMINA', 'snpeff', 'snpEff -version', containerDir.resolve('snpeff:5.0.sif'))
    }
    if (params.virus == 'sars-cov2') {
        specs << toolSpec('ILLUMINA', 'pangolin', 'pangolin --version', containerDir.resolve('pangolin:4.4.sif'))
        specs << toolSpec('ILLUMINA', 'nextclade', 'nextclade --version', containerDir.resolve('nextclade:3.18.sif'))
    }
    specs
}

def metadataDir(params) {
    java.nio.file.Path.of(params.outDir.toString()).toAbsolutePath().normalize()
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
    catch (Exception ignored) {
        return value.toString()
    }
}

def safeMetadataValue(closure) {
    try {
        return closure.call()
    }
    catch (Exception ignored) {
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
    catch (Exception ignored) {
        return null
    }
}

def metadataFailureMessage(workflow) {
    safeMetadataValue { -> workflow.errorMessage } ?: safeMetadataValue { -> workflow.errorReport }
}

def writeRunManifest(workflow, params, status, failureMessage = null) {
    def outputDir = metadataDir(params)
    java.nio.file.Files.createDirectories(outputDir)

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
            input_dir: absoluteMetadataPath(params.inDir),
            output_dir: absoluteMetadataPath(params.outDir)
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
            software_versions: 'software_versions.tsv',
            containers: 'container_manifest.tsv',
            trace: 'execution_trace.tsv',
            report: 'execution_report.html',
            timeline: 'execution_timeline.html'
        ]
    ]

    def target = outputDir.resolve('run_manifest.json')
    def temporary = outputDir.resolve('run_manifest.json.tmp')
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

def containerSpecs(params, workflow) {
    def specs = []
    if (params.mode == 'NANOPORE') {
        specs << localContainerSpec('nanopore_base', params.base_container)
        specs << remoteContainerSpec('clair3', 'docker://hkubal/clair3:v1.2.0')
    }
    else if (params.mode == 'ILLUMINA') {
        def containerDir = java.nio.file.Path.of(workflow.projectDir.toString()).resolve('containers')
        specs.addAll([
            localContainerSpec('edirect', containerDir.resolve('edirect:1.1.0.sif')),
            localContainerSpec('generate_consensus', containerDir.resolve('generate_consensus:2.0.0.sif')),
            localContainerSpec('fastp', containerDir.resolve('fastp:1.0.1.sif')),
            localContainerSpec('samtools', containerDir.resolve('samtools:1.11.0.sif')),
            localContainerSpec('mafft', containerDir.resolve('mafft:7.505_2.sif')),
            localContainerSpec('picard', containerDir.resolve('picard:2.27.2_2.sif')),
            localContainerSpec('intrahost_analysis', containerDir.resolve('intrahost_analysis:1.1.0.sif')),
            localContainerSpec('generate_plots', containerDir.resolve('generate_plots:2.0.0.sif')),
            localContainerSpec('compiled_outputs', containerDir.resolve('compiled_outputs:1.1.0.sif'))
        ])
        if (params.runSnpEff) {
            specs << localContainerSpec('snpeff', containerDir.resolve('snpeff:5.0.sif'))
            specs << localContainerSpec('generate_report', containerDir.resolve('generate_report:1.1.0.sif'))
        }
        if (params.virus == 'sars-cov2') {
            specs << localContainerSpec('pangolin', containerDir.resolve('pangolin:4.4.sif'))
            specs << localContainerSpec('nextclade', containerDir.resolve('nextclade:3.18.sif'))
        }
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
            toolSpec('NANOPORE', 'clair3', 'run_clair3.sh -v ', 'docker://hkubal/clair3:v1.2.0')
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

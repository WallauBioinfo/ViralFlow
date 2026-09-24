include {normalizeTrimLen} from './param_helpers.nf'

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

// Where Nextflow writes one of its own reports - `trace`, `report` or
// `timeline` - as the session resolved it, or null when that report is off.
//
// Read from the session rather than derived from params.outDir, because the two
// can disagree. nextflow.config builds these paths from params.outDir while it
// is parsed, before a -c config is merged, so an outDir set there moves every
// other output but not these three. -with-report, -with-timeline and
// -with-trace move them anywhere; nf-test always does so for the trace. A
// manifest that assumed outDir would name files that are not there.
//
// A relative path is resolved against the launch directory, as Nextflow does.
def nextflowReportPath(workflow, scope) {
    safeMetadataValue { ->
        def options = nextflow.Global.session.config[scope] as Map
        options?.enabled && options.file
            ? absoluteMetadataPath(workflow.launchDir.resolve(options.file.toString()))
            : null
    }
}

// The executor Nextflow runs a process on when the process does not choose its
// own, resolved the way Nextflow resolves it (ExecutorFactory.getExecutorName,
// as of 26.04.6): process.executor, then executor.name, then the NXF_EXECUTOR
// environment variable, then local. -process.executor=... on the command line
// sets the first. A withName or withLabel selector can still send one process
// elsewhere; no ViralFlow configuration does.
//
// Read from the session config rather than guessed from the profile name,
// which recorded any executor but PBS as local.
def configuredExecutor() {
    safeMetadataValue { ->
        def config = nextflow.Global.session.config
        def processExecutor = (config.process as Map)?.executor
        def executorName = (config.executor as Map)?.name
        (processExecutor ?: (executorName instanceof String ? executorName : null)
            ?: System.getenv('NXF_EXECUTOR') ?: 'local').toString()
    }
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
            executor: configuredExecutor(),
            // The engine Nextflow actually ran tasks with, not a guess from the
            // profile name: that recorded every -profile docker run as
            // singularity. Null when no container engine is enabled.
            container_engine: safeMetadataValue { -> workflow.containerEngine }
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
            trace: nextflowReportPath(workflow, 'trace'),
            report: nextflowReportPath(workflow, 'report'),
            timeline: nextflowReportPath(workflow, 'timeline')
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

// The reference a container engine accepts for an image. Singularity and
// Apptainer need the docker:// prefix on a registry image; the Docker engine
// rejects it ("invalid reference format"). Strip it for Docker only - the rule
// the runClair3 directive in nextflow.config applies, which a config file
// cannot share by calling this function.
def engineImageReference(reference, engine) {
    def value = reference.toString()
    engine == 'docker' ? value.replaceFirst(/^docker:\/\//, '') : value
}

// An image the Docker daemon resolves: a tag or digest reference, never a file.
// base_container is one under -profile docker (viralflow/nanopore-base:...),
// and recording it as a local SIF turned the tag into a path that does not exist.
def dockerImageSpec(name, reference) {
    [name: name, kind: 'docker_image', identity: engineImageReference(reference, 'docker')]
}

def toolSpec(mode, name, command, containerValue) {
    [mode: mode, tool: name, command: command, container: containerValue.toString()]
}

// The channel builders below own the spec-map -> tuple mapping. captureToolVersion
// and captureContainerMetadata read those tuples positionally, so keeping the
// mapping in one place is what lets a test pin the field order.
//
// engine is the container engine actually in use: main.nf passes
// workflow.containerEngine. An argument rather than read here so a test can pin
// it; the suite runs under -profile docker in CI and singularity elsewhere.
def containerSpecChannel(params, engine) {
    channel.fromList(
        containerSpecs(params, engine).collect { spec ->
            tuple(spec.name, spec.kind, spec.identity)
        }
    )
}

def toolSpecChannel(params, workflow, engine) {
    channel.fromList(
        toolSpecs(params, workflow, engine).collect { spec ->
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

def containerSpecs(params, engine) {
    def specs = []
    if (params.mode == 'NANOPORE') {
        if (engine == 'docker') {
            // Both are references the Docker daemon resolves, not files: the
            // docker profile sets base_container to an image tag, and Clair3 is
            // pulled by digest.
            specs << dockerImageSpec('nanopore_base', params.base_container)
            specs << dockerImageSpec('clair3', params.clair3_container)
        }
        else {
            specs << localContainerSpec('nanopore_base', params.base_container)
            specs << remoteContainerSpec('clair3', params.clair3_container)
        }
        // These two are the whole list: GENPLOTS runs after NANOPORE too, but
        // in NANOPORE mode its processes use the base image (see
        // configs/containers.config). An ILLUMINA image listed here would be
        // one the run never started.
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

// Each version command runs in the image it reports on, so the container must be
// a reference the engine in use accepts; see engineImageReference().
def toolSpecs(params, workflow, engine) {
    if (params.mode == 'NANOPORE') {
        def base = engineImageReference(params.base_container, engine)
        def clair3 = engineImageReference(params.clair3_container, engine)
        def specs = [
            toolSpec('NANOPORE', 'porechop_abi', 'porechop_abi --version', base),
            toolSpec('NANOPORE', 'minimap2', 'minimap2 --version', base),
            toolSpec('NANOPORE', 'samtools', 'samtools --version | head -n 1', base),
            toolSpec('NANOPORE', 'bcftools', 'bcftools --version | head -n 1', base),
            toolSpec('NANOPORE', 'clair3', 'run_clair3.sh -v ', clair3),
            // Draws the GENPLOTS coverage plot, which every run executes in the
            // base image. GENPLOTS' samtools is the one already listed above.
            toolSpec('NANOPORE', 'bamdash', 'bamdash --version', base)
        ]
        // Only when --trimLen asks for it, from the same rule NANOPORE.nf uses
        // to decide whether runBamUtils runs at all. Listing it unconditionally
        // would report a tool that never touched the reads; omitting it when
        // trimming is on leaves the step that rewrote every alignment out of the
        // provenance record.
        //
        // bamUtil has no --version flag: `bam` alone exits 255 without printing
        // one, and `bam help` puts it on the second line. sed rather than head
        // because `bam help` keeps writing afterwards, and head closing the pipe
        // early makes the command exit 141 under pipefail.
        if (normalizeTrimLen(params.trimLen) > 0) {
            specs << toolSpec('NANOPORE', 'bamutil', 'bam help 2>&1 | sed -n 2p', base)
        }
        // samtools already covers --primersBED: ampliconclip is a samtools
        // subcommand, so that option adds no tool of its own.
        return specs
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

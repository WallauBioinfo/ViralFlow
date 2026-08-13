import groovy.json.JsonOutput

import java.nio.file.Files
import java.nio.file.Path

class MetadataHelper {

    static Path metadataDir(def params) {
        return Path.of(params.outDir.toString()).toAbsolutePath().normalize()
            .resolve("RUN_METADATA")
    }

    static void writeManifest(def workflow, def params, String status, String failureMessage = null) {
        Path outputDir = metadataDir(params)
        Files.createDirectories(outputDir)

        def manifest = [
            schema_version: 1,
            pipeline      : [
                name      : safeValue { workflow.manifest.name },
                version   : safeValue { workflow.manifest.version },
                repository: safeValue { workflow.manifest.homePage },
                revision  : gitValue(workflow.projectDir, ["rev-parse", "--abbrev-ref", "HEAD"]),
                commit_id : gitValue(workflow.projectDir, ["rev-parse", "HEAD"]),
                git_dirty : gitDirty(workflow.projectDir)
            ],
            execution     : [
                status          : status,
                failure_message : failureMessage,
                session_id      : safeValue { workflow.sessionId },
                run_name        : safeValue { workflow.runName },
                command_line    : safeValue { workflow.commandLine },
                profile         : safeValue { workflow.profile },
                nextflow_version: nextflow.BuildInfo.version,
                start_time      : normalize(safeValue { workflow.start }),
                end_time        : status == "RUNNING" ? null : normalize(new Date()),
                duration        : normalize(safeValue { workflow.duration }),
                success         : status == "SUCCESS"
            ],
            runtime       : [
                user        : safeValue { workflow.userName },
                host        : hostname(),
                os          : System.getProperty("os.name"),
                os_version  : System.getProperty("os.version"),
                architecture: System.getProperty("os.arch"),
                executor    : executorName(safeValue { workflow.profile }),
                container_engine: containerEngine(safeValue { workflow.profile })
            ],
            paths         : [
                launch_dir : absolute(safeValue { workflow.launchDir }),
                project_dir: absolute(safeValue { workflow.projectDir }),
                work_dir   : absolute(safeValue { workflow.workDir }),
                input_dir  : absolute(params.inDir),
                output_dir : absolute(params.outDir)
            ],
            analysis      : [
                mode        : normalize(params.mode),
                virus       : normalize(params.virus),
                clair3_model: params.mode == "NANOPORE" ? normalize(params.clair3_model) : null,
                clair3_qual : params.mode == "NANOPORE" ? normalize(params.clair3_qual) : null,
                mapping_quality: params.mode == "NANOPORE" ? normalize(params.mapping_quality) : null,
                af_threshold: params.mode == "NANOPORE" ? normalize(params.af_threshold) : null,
                min_depth   : params.mode == "NANOPORE" ? normalize(params.np_min_depth) : normalize(params.depth),
                consensus_mask_rule: params.mode == "NANOPORE" ? "depth <= min_depth" : null
            ],
            parameters    : normalize(params.entrySet().collectEntries { entry ->
                [(entry.key.toString()): entry.value]
            }),
            files         : [
                input_checksums : "input_checksums.tsv",
                software_versions: "software_versions.tsv",
                containers      : "container_manifest.tsv",
                trace           : "execution_trace.tsv",
                report          : "execution_report.html",
                timeline        : "execution_timeline.html"
            ]
        ]

        Path target = outputDir.resolve("run_manifest.json")
        Path temporary = outputDir.resolve("run_manifest.json.tmp")
        temporary.toFile().text = JsonOutput.prettyPrint(JsonOutput.toJson(manifest)) + System.lineSeparator()
        Files.move(
            temporary,
            target,
            java.nio.file.StandardCopyOption.REPLACE_EXISTING,
            java.nio.file.StandardCopyOption.ATOMIC_MOVE
        )
    }

    static List<Map> containerSpecs(def params, def workflow) {
        def specs = []

        if (params.mode == "NANOPORE") {
            specs << localContainerSpec("nanopore_base", params.base_container)
            specs << remoteContainer("clair3", "docker://hkubal/clair3:v1.2.0")
        } else if (params.mode == "ILLUMINA") {
            def projectDir = Path.of(workflow.projectDir.toString())
            specs.addAll([
                localContainerSpec("edirect", projectDir.resolve("containers/edirect:1.1.0.sif")),
                localContainerSpec("generate_consensus", projectDir.resolve("containers/generate_consensus:2.0.0.sif")),
                localContainerSpec("fastp", projectDir.resolve("containers/fastp:1.0.1.sif")),
                localContainerSpec("samtools", projectDir.resolve("containers/samtools:1.11.0.sif")),
                localContainerSpec("mafft", projectDir.resolve("containers/mafft:7.505_2.sif")),
                localContainerSpec("picard", projectDir.resolve("containers/picard:2.27.2_2.sif")),
                localContainerSpec("intrahost_analysis", projectDir.resolve("containers/intrahost_analysis:1.1.0.sif")),
                localContainerSpec("generate_plots", projectDir.resolve("containers/generate_plots:2.0.0.sif")),
                localContainerSpec("compiled_outputs", projectDir.resolve("containers/compiled_outputs:1.1.0.sif"))
            ])
            if (params.runSnpEff) {
                specs << localContainerSpec("snpeff", projectDir.resolve("containers/snpeff:5.0.sif"))
                specs << localContainerSpec("generate_report", projectDir.resolve("containers/generate_report:1.1.0.sif"))
            }
            if (params.virus == "sars-cov2") {
                specs << localContainerSpec("pangolin", projectDir.resolve("containers/pangolin:4.4.sif"))
                specs << localContainerSpec("nextclade", projectDir.resolve("containers/nextclade:3.18.sif"))
            }
        }

        return specs.unique { it.identity }
    }

    static List<Map> toolSpecs(def params, def workflow) {
        if (params.mode == "NANOPORE") {
            return [
                tool("NANOPORE", "porechop_abi", "porechop_abi --version", params.base_container),
                tool("NANOPORE", "minimap2", "minimap2 --version", params.base_container),
                tool("NANOPORE", "samtools", "samtools --version | head -n 1", params.base_container),
                tool("NANOPORE", "bcftools", "bcftools --version | head -n 1", params.base_container),
                tool("NANOPORE", "clair3", "run_clair3.sh -v ", "docker://hkubal/clair3:v1.2.0")
            ]
        }

        def projectDir = Path.of(workflow.projectDir.toString()).resolve("containers")
        def specs = [
            tool("ILLUMINA", "fastp", "fastp --version", projectDir.resolve("fastp:1.0.1.sif")),
            tool("ILLUMINA", "bwa", "bwa 2>&1 | head -n 3", projectDir.resolve("generate_consensus:2.0.0.sif")),
            tool("ILLUMINA", "samtools", "samtools --version | head -n 1", projectDir.resolve("generate_consensus:2.0.0.sif")),
            tool("ILLUMINA", "ivar", "ivar version", projectDir.resolve("generate_consensus:2.0.0.sif")),
            tool("ILLUMINA", "mafft", "mafft --version", projectDir.resolve("mafft:7.505_2.sif"))
        ]
        if (params.runSnpEff) {
            specs << tool("ILLUMINA", "snpeff", "snpEff -version", projectDir.resolve("snpeff:5.0.sif"))
        }
        if (params.virus == "sars-cov2") {
            specs << tool("ILLUMINA", "pangolin", "pangolin --version", projectDir.resolve("pangolin:4.4.sif"))
            specs << tool("ILLUMINA", "nextclade", "nextclade --version", projectDir.resolve("nextclade:3.18.sif"))
        }
        return specs
    }

    static String failureMessage(def workflow) {
        return safeValue { workflow.errorMessage } ?: safeValue { workflow.errorReport }
    }

    static Map normalizeMap(Map value) {
        return normalize(value) as Map
    }

    static Map localContainerSpec(String name, def pathValue) {
        String identity = absolute(pathValue)
        String kind = Files.isDirectory(Path.of(identity)) ? "local_sandbox" : "local_sif"
        return [name: name, kind: kind, identity: identity]
    }

    private static Map remoteContainer(String name, String uri) {
        return [name: name, kind: "remote_uri", identity: uri]
    }

    private static Map tool(String mode, String name, String command, def containerValue) {
        return [
            mode     : mode,
            tool     : name,
            command  : command,
            container: containerValue.toString()
        ]
    }

    private static Object normalize(def value) {
        if (value == null || value instanceof Boolean || value instanceof Number || value instanceof String) {
            return value
        }
        if (value instanceof Path || value instanceof File) {
            return absolute(value)
        }
        if (value instanceof Date) {
            return value.toInstant().toString()
        }
        if (value instanceof Map) {
            return value.collectEntries { key, item -> [(key.toString()): normalize(item)] }
        }
        if (value instanceof Collection) {
            return value.collect { normalize(it) }
        }
        if (value.getClass().isArray()) {
            return value.toList().collect { normalize(it) }
        }
        return value.toString()
    }

    private static String absolute(def value) {
        if (value == null) {
            return null
        }
        try {
            return Path.of(value.toString()).toAbsolutePath().normalize().toString()
        } catch (Exception ignored) {
            return value.toString()
        }
    }

    private static def safeValue(Closure closure) {
        try {
            return closure.call()
        } catch (Exception ignored) {
            return null
        }
    }

    private static String gitValue(def projectDir, List<String> arguments) {
        try {
            def command = ["git", "-C", projectDir.toString()] + arguments
            def process = new ProcessBuilder(command).redirectErrorStream(true).start()
            String output = process.inputStream.text.trim()
            return process.waitFor() == 0 ? output : null
        } catch (Exception ignored) {
            return null
        }
    }

    private static Boolean gitDirty(def projectDir) {
        String status = gitValue(projectDir, ["status", "--porcelain"])
        return status == null ? null : !status.isEmpty()
    }

    private static String hostname() {
        try {
            return InetAddress.localHost.hostName
        } catch (Exception ignored) {
            return null
        }
    }

    private static String containerEngine(def profileValue) {
        String profile = profileValue?.toString() ?: ""
        if (profile.contains("apptainer")) {
            return "apptainer"
        }
        return "singularity"
    }

    private static String executorName(def profileValue) {
        String profile = profileValue?.toString() ?: ""
        return profile.contains("pbs") ? "pbs" : "local"
    }
}

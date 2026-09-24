nextflow.enable.dsl = 2

include {
    METADATA
    captureContainerMetadata
    captureContainerMetadata as capture_missing_container
} from '../../modules/metadata.nf'
include {
    localContainerSpec
    normalizeMetadata
    referenceMetadataChannel
    containerSpecChannel
    toolSpecChannel
} from '../../modules/metadata_helpers.nf'
include { processInputs } from '../../workflows/step0-input-handling.nf'

// The container engine the spec helpers classify for. Pinned by a test's params
// rather than read from workflow.containerEngine, so a test asserts the same
// thing under -profile docker in CI as under singularity elsewhere. main.nf
// passes the real engine.
def fixtureEngine() {
    params.getOrDefault('engine', 'singularity')
}

workflow METADATA_FIXTURE {
    main:
        def inputFile = file("${projectDir}/tests/data/bcftools/ref.fa")
        // Take the container from the parameter rather than hardcoding the SIF
        // path, so this runs under whichever profile is active - the .sif under
        // singularity, the image under -profile docker.
        def containerPath = params.base_container

        checksum_inputs = channel.of(
            tuple(
                "reference",
                "reference_fasta",
                inputFile.toAbsolutePath().normalize().toString(),
                inputFile
            )
        )
        tool_specs = channel.of(
            tuple(
                "NANOPORE",
                "bcftools",
                "bcftools --version | head -n 1",
                containerPath.toString()
            )
        )
        container_specs = channel.of(
            tuple("test_remote", "remote_uri", "docker://example/test:1.0"),
            tuple(
                "test_sandbox",
                "local_sandbox",
                file("${projectDir}/tests/data/metadata/sandbox-container")
                    .toAbsolutePath()
                    .normalize()
                    .toString()
            ),
            tuple(
                "test_sif",
                "local_sif",
                inputFile.toAbsolutePath().normalize().toString()
            )
        )

        METADATA(
            checksum_inputs,
            tool_specs,
            container_specs,
            channel.of([["sample", 1, "single", inputFile.toString(), ""]]),
            params.metadataDir
        )

    emit:
        checksums = METADATA.out.input_checksums
        versions = METADATA.out.software_versions
        containers = METADATA.out.container_manifest
}

workflow MISSING_CONTAINER_FIXTURE {
    main:
        missing_spec = localContainerSpec(
            "missing",
            file("${projectDir}/tests/data/metadata/missing.sif")
        )
        missingContainerCh = channel.of(
            tuple(
                missing_spec.name,
                missing_spec.kind,
                missing_spec.identity
            )
        )

        capture_missing_container(missingContainerCh)
}

workflow CLASSIFY_LOCAL_CONTAINERS_FIXTURE {
    main:
        specs = [
            localContainerSpec(
                "sandbox",
                file("${projectDir}/tests/data/metadata/sandbox-container")
            ),
            localContainerSpec(
                "sif",
                file("${projectDir}/tests/data/bcftools/ref.fa")
            )
        ]
        specsCh = channel.value(specs)

    emit:
        specsCh
}

workflow NORMALIZE_METADATA_FIXTURE {
    main:
        normalized = normalizeMetadata([
            null_value: null,
            enabled: true,
            count: 3,
            memory: 4.GB,
            values: ["a", 2],
            path: file("${projectDir}/tests/data/bcftools/ref.fa")
        ])
        normalizedCh = channel.value(normalized)

    emit:
        normalizedCh
}

// containerSpecs() and toolSpecs() build the tuples that captureContainerMetadata
// and captureToolVersion consume by position. Neither builder was executed by any
// test: main.nf reaches them only after processInputs(), and the one test that runs
// main.nf aborts in validation first. A field reorder would therefore corrupt
// container_manifest.tsv with every test still green.
//
// This drives the real builder into the real process, so the two stay in agreement.
workflow CONTAINER_SPECS_FIXTURE {
    main:
        captureContainerMetadata(containerSpecChannel(params, fixtureEngine()))

    emit:
        rows = captureContainerMetadata.out
}

// The ILLUMINA branch of containerSpecs() cannot go through
// captureContainerMetadata the way CONTAINER_SPECS_FIXTURE does: that process
// checksums each local .sif, and the ILLUMINA images are pulled by
// `viralflow build-containers` rather than living in the repository, so they are
// absent wherever the suite runs. Emit the tuples instead - enough to pin which
// containers the branch declares and how each is classified.
workflow CONTAINER_SPECS_TUPLES_FIXTURE {
    main:
        containerSpecsCh = containerSpecChannel(params, fixtureEngine())

    emit:
        containerSpecsCh
}

// toolSpecs feeds captureToolVersion, which can only run inside each tool's own
// container. Assert the tuple contract here; the execution path is covered by
// METADATA_FIXTURE.
workflow TOOL_SPECS_FIXTURE {
    main:
        toolSpecsCh = toolSpecChannel(params, workflow, fixtureEngine())

    emit:
        toolSpecsCh
}

// Drives the real processInputs outputs through the same helper main.nf uses,
// so a regression in reference metadata handling fails here instead of hiding
// behind a synthetic channel.
workflow OPTIONAL_GFF_METADATA_FIXTURE {
    main:
        processInputs()

        fastaMetadataCh = referenceMetadataChannel(
            "reference_fasta",
            processInputs.out.refFa
        )
        gffMetadataCh = params.mode == "ILLUMINA"
            ? referenceMetadataChannel(
                "referenceGff",
                processInputs.out.refGff
            )
            : channel.empty()

    emit:
        gffMetadataCh
        fastaMetadataCh
}

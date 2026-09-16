nextflow.enable.dsl = 2

include {
    METADATA
    capture_container_metadata as capture_missing_container
} from '../../modules/metadata.nf'
include {
    localContainerSpec
    normalizeMetadata
    referenceMetadataChannel
} from '../../modules/metadata_helpers.nf'
include { processInputs } from '../../workflows/step0-input-handling.nf'

workflow METADATA_FIXTURE {
    main:
        def inputFile = file("${projectDir}/tests/data/bcftools/ref.fa")
        def containerPath = file("${projectDir}/containers/baseContainer.sif")

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
        missing_container_ch = channel.of(
            tuple(
                missing_spec.name,
                missing_spec.kind,
                missing_spec.identity
            )
        )

        capture_missing_container(missing_container_ch)
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
        specs_ch = channel.value(specs)

    emit:
        specs_ch
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
        normalized_ch = channel.value(normalized)

    emit:
        normalized_ch
}

// Drives the real processInputs outputs through the same helper main.nf uses,
// so a regression in reference metadata handling fails here instead of hiding
// behind a synthetic channel.
workflow OPTIONAL_GFF_METADATA_FIXTURE {
    main:
        processInputs()

        fasta_metadata_ch = referenceMetadataChannel(
            "reference_fasta",
            processInputs.out.ref_fa
        )
        gff_metadata_ch = params.mode == "ILLUMINA"
            ? referenceMetadataChannel(
                "reference_gff",
                processInputs.out.ref_gff
            )
            : channel.empty()

    emit:
        gff_metadata_ch
        fasta_metadata_ch
}

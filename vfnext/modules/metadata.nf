process checksum_metadata_input {
    tag "${sample_id}:${role}"

    input:
        tuple val(sample_id), val(role), val(original_path), path(input_file)

    output:
        path("${sample_id}.${role}.checksum.tsv")

    script:
    """
    set -euo pipefail

    checksum=\$(sha256sum ${input_file} | cut -d ' ' -f 1)
    size=\$(stat -Lc '%s' ${input_file})
    printf '%s\\t%s\\t%s\\t%s\\t%s\\n' \
        '${sample_id}' '${role}' '${original_path}' "\${size}" "\${checksum}" \
        > ${sample_id}.${role}.checksum.tsv
    """
}

process capture_tool_version {
    tag "${tool_name}"
    container "${container_identity}"

    input:
        tuple val(mode), val(tool_name), val(version_command), val(container_identity)

    output:
        path("${tool_name}.version.tsv")

    script:
    """
    set -euo pipefail

    set +e
    raw=\$(bash -o pipefail -c '${version_command}' 2>&1)
    status=\$?
    set -e

    cleaned=\$(printf '%s' "\${raw}" | tr '\\t\\r\\n' '   ' | sed 's/  */ /g; s/^ //; s/ \$//')
    if [[ -z "\${cleaned}" ]]; then
        echo "Unable to determine ${tool_name} version" >&2
        exit 1
    fi

    version=\$(printf '%s' "\${cleaned}" | sed -E 's/.*([vV]?[0-9]+([.][0-9A-Za-z_-]+)+).*/\\1/' | awk '{print \$1}')
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
        '${mode}' '${tool_name}' "\${version}" '${container_identity}' "\${status}" "\${cleaned}" \
        > ${tool_name}.version.tsv
    """
}

process capture_container_metadata {
    tag "${container_name}"

    input:
        tuple val(container_name), val(container_kind), val(container_identity)

    output:
        path("${container_name}.container.tsv")

    script:
    """
    set -euo pipefail

    case '${container_kind}' in
        local_sif)
        if [[ -d '${container_identity}' ]]; then
            echo "Configured SIF path is a directory: ${container_identity}" >&2
            exit 1
        fi
        if [[ ! -f '${container_identity}' ]]; then
            echo "Configured SIF file does not exist: ${container_identity}" >&2
            exit 1
        fi
        checksum=\$(sha256sum '${container_identity}' | cut -d ' ' -f 1)
        size=\$(stat -Lc '%s' '${container_identity}')
            ;;
        local_sandbox)
            if [[ ! -d '${container_identity}' ]]; then
                echo "Configured sandbox directory does not exist: ${container_identity}" >&2
                exit 1
            fi
        checksum='NA'
        size='NA'
            ;;
        remote_uri)
            checksum='NA'
            size='NA'
            ;;
        *)
            echo "Unsupported container metadata kind: ${container_kind}" >&2
            exit 1
            ;;
    esac

    printf '%s\\t%s\\t%s\\t%s\\t%s\\n' \
        '${container_name}' '${container_kind}' '${container_identity}' "\${size}" "\${checksum}" \
        > ${container_name}.container.tsv
    """
}

workflow METADATA {
    take:
        checksum_inputs
        tool_specs
        container_specs
        resolved_inputs
        metadata_dir

    main:
        checksum_metadata_input(checksum_inputs)
        capture_tool_version(tool_specs)
        capture_container_metadata(container_specs)

        resolved_inputs
            .flatMap { rows -> rows }
            .map { row -> row.collect { value -> value.toString().replace('\t', ' ').replace('\n', ' ') }.join('\t') + '\n' }
            .collectFile(
                name: "resolved_sample_inputs.tsv",
                seed: "sample_id\\tchunk_index\\tlayout\\tfastq_1\\tfastq_2\\n",
                sort: false,
                newLine: false,
                storeDir: metadata_dir
            )
            .set { resolved_sample_inputs }

        checksum_metadata_input.out
            .collectFile(
                name: "input_checksums.tsv",
                seed: "sample_id\trole\tabsolute_path\tsize_bytes\tsha256\n",
                sort: true,
                newLine: false,
                storeDir: metadata_dir
            )
            .set { input_checksums }

        capture_tool_version.out
            .collectFile(
                name: "software_versions.tsv",
                seed: "mode\ttool\tversion\tcontainer\texit_status\traw_output\n",
                sort: true,
                newLine: false,
                storeDir: metadata_dir
            )
            .set { software_versions }

        capture_container_metadata.out
            .collectFile(
                name: "container_manifest.tsv",
                seed: "name\tkind\tidentity\tsize_bytes\tsha256\n",
                sort: true,
                newLine: false,
                storeDir: metadata_dir
            )
            .set { container_manifest }

    emit:
        input_checksums
        software_versions
        container_manifest
        resolved_sample_inputs
}

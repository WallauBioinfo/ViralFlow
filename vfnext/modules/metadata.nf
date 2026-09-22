// Shared shell helpers for the two metadata processes that run on the host
// rather than inside a container. GNU coreutils is not guaranteed there:
// `sha256sum` and `stat -Lc` do not exist on macOS/BSD, where the original
// commands aborted the whole run before any analysis started.
//
// These processes cannot simply be containerised. captureContainerMetadata
// receives the container identity as a value, not a staged path, so the file
// it must checksum is not visible from inside a container.
//
// `wc -c < file` is POSIX and follows symlinks via the shell redirect, giving
// the same result as `stat -L`. BSD `wc` pads its output, hence the trim.
def portableFileMetrics() {
    return '''
sha256_checksum() {
    if command -v sha256sum > /dev/null 2>&1; then
        sha256sum "$1" | cut -d ' ' -f 1
    elif command -v shasum > /dev/null 2>&1; then
        shasum -a 256 "$1" | cut -d ' ' -f 1
    else
        echo "No SHA-256 utility found; need sha256sum or shasum on PATH" >&2
        return 1
    fi
}

file_size_bytes() {
    wc -c < "$1" | tr -d '[:space:]'
}
'''
}

process checksumMetadataInput {
    tag "${sample_id}:${role}"

    input:
        tuple val(sample_id), val(role), val(original_path), path(input_file)

    output:
        path("${sample_id}.${role}.checksum.tsv")

    script:
    """
    set -euo pipefail
    ${portableFileMetrics()}
    checksum=\$(sha256_checksum ${input_file})
    size=\$(file_size_bytes ${input_file})
    printf '%s\\t%s\\t%s\\t%s\\t%s\\n' \
        '${sample_id}' '${role}' '${original_path}' "\${size}" "\${checksum}" \
        > ${sample_id}.${role}.checksum.tsv
    """
}

process captureToolVersion {
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

    # Leftmost match, via bash's own =~, rather than a greedy sed expression.
    if [[ "\${cleaned}" =~ ([vV]?[0-9]+([.][0-9A-Za-z_-]+)+) ]]; then
        version="\${BASH_REMATCH[1]}"
    else
        # Nothing version-shaped in the output: keep the cleaned text
        version="\${cleaned}"
    fi

    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
        '${mode}' '${tool_name}' "\${version}" '${container_identity}' "\${status}" "\${cleaned}" \
        > ${tool_name}.version.tsv
    """
}

process captureContainerMetadata {
    tag "${container_name}"

    input:
        tuple val(container_name), val(container_kind), val(container_identity)

    output:
        path("${container_name}.container.tsv")

    script:
    """
    set -euo pipefail
    ${portableFileMetrics()}
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
            checksum=\$(sha256_checksum '${container_identity}')
            size=\$(file_size_bytes '${container_identity}')
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
        resolvedInputs
        metadata_dir

    main:
        checksumMetadataInput(checksum_inputs)
        captureToolVersion(tool_specs)
        captureContainerMetadata(container_specs)

        resolvedInputs
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

        checksumMetadataInput.out
            .collectFile(
                name: "input_checksums.tsv",
                seed: "sample_id\trole\tabsolute_path\tsize_bytes\tsha256\n",
                sort: true,
                newLine: false,
                storeDir: metadata_dir
            )
            .set { input_checksums }

        captureToolVersion.out
            .collectFile(
                name: "software_versions.tsv",
                seed: "mode\ttool\tversion\tcontainer\texit_status\traw_output\n",
                sort: true,
                newLine: false,
                storeDir: metadata_dir
            )
            .set { software_versions }

        captureContainerMetadata.out
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

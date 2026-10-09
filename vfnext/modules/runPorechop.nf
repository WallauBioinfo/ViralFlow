
process runPorechop {
    label "NP_basecontainer"
    // Define the process parameters
    publishDir { "${params.outDir}/${meta.id}_results/" }, mode: 'copy', overwrite: true, pattern: "*.fastq.gz"
    tag "${meta.id}"

    input:
        tuple val(meta), path(fastq)

    output:
        tuple val(meta), path("${meta.id}.chopped.fastq.gz")

    script:
    // Compressed before publishing: uncompressed, the trimmed reads are some
    // 35 times larger (6.9 MB against 188 KB on the truth fixture). minimap2
    // reads the BGZF output directly.
    //
    // Porechop writes a plain file that bgzip then replaces, rather than
    // streaming its stdout into bgzip: -abi runs a helper that prints its
    // progress to stdout ahead of the reads, which would corrupt the FASTQ.
    // Porechop's own .gz output also goes through an uncompressed temporary
    // file, and compresses it with single-threaded gzip.
    """
    set -euo pipefail

    porechop_abi -abi -i ${fastq} -t ${task.cpus} -o ${meta.id}.chopped.fastq
    bgzip -@ ${task.cpus} ${meta.id}.chopped.fastq
    """
}

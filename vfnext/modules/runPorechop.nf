
process run_porechop {
    label "NP_basecontainer"
    // Define the process parameters
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true
    tag "${meta.id}"

    input:
        tuple val(meta), path(fastq)
    
    output:
        tuple val(meta), path("${meta.id}.chopped.fastq")
  
    script:
    """
    porechop_abi -abi -i ${fastq} -t ${task.cpus} -o ${meta.id}.chopped.fastq
    """
}
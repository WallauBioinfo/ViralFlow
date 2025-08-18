process run_minimap2 {
    label "NP_basecontainer"
    // Define the process parameters
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true
    tag "${meta.id}"

    input:
        tuple val(meta), path(fastq)
        path(ref)
    
    output:
        tuple val(meta), path("${meta.id}.sorted.bam*")
  
    script:
    """
    minimap2 -a -x map-ont -t ${task.cpus} ${ref} ${fastq} \
    | samtools view -bS -F 4 - \
    | samtools sort -o ${meta.id}.sorted.bam
        
    samtools index ${meta.id}.sorted.bam
    """
}

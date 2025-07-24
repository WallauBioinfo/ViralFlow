
process run_bam_utils {
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true
    tag "${meta.id}"

    input: 
        tuple val(meta), path(bam)
    
    output:
        tuple val(meta), path("${meta.id}.trim.sorted.bam*")

    script:
    """
    bam trimBam ${bam[0]} ${meta.id}.trim.bam -L ${params.trim_len} -R ${params.trim_len}
    samtools sort ${meta.id}.trim.bam -o ${meta.id}.trim.sorted.bam
    samtools index ${meta.id}.trim.sorted.bam
    """
}

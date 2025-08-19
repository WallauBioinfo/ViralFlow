
process run_amplicon_clip {
    // Define the process parameters
    publishDir "${params.outDir}/${meta.id}", mode: 'copy', overwrite: true
    tag "${meta.id}"

    input:
        tuple val(meta), path(sorted_bam)
        path(primer_bed)
    
    output:
        tuple val(meta), path("${meta.id}.primer_clip.bam*")
  
    script:
    """
    samtools ampliconclip --strand --hard-clip -b ${primer_bed} ${sorted_bam} -f ./trimmed_reads.txt | samtools sort -o ${meta.id}.primer_clip.bam
    samtools index ${meta.id}.primer_clip.bam
    """
}

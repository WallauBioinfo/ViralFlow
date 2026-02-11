
process align2ref{
  publishDir "${params.outDir}/${meta.id}_results/", mode: "copy"
  
  input:
    tuple val(meta), path(reads), val(is_paired_end), path(fasta_amb), path(fasta_ann), path(fasta_bwt), path(fasta_pac), path(fasta_sa)
    path(ref_fa)

  output:
    tuple val(meta), path("*.sorted.bam"), path("*.bai"), val(is_paired_end), emit: regular_output
    
  script:
    sample_id = meta.id
    """
    # Link reference files
    if [[ ! -f ./${fasta_amb.getSimpleName()}.fasta ]]; then
        ln -s ${ref_fa} ./${fasta_amb.getSimpleName()}.fasta
    fi

    if [[ ${is_paired_end} == true ]]; then
        bwa mem ./${ref_fa} ${reads[0]} ${reads[1]} \
                -o ${sample_id}.bam -t ${params.bwa_threads} 
    else
        bwa mem ./${ref_fa} ${reads[0]} \
                -o ${sample_id}.bam -t ${params.bwa_threads} 
    fi

    # Sort and index
    samtools sort -o ${sample_id}.sorted.bam ${sample_id}.bam
    samtools index ${sample_id}.sorted.bam
    """
}

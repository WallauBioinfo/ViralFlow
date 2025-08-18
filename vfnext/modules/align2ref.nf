
process align2ref{
  tag "${meta.id}"
  publishDir "${params.outDir}/${meta.id}_results/", mode: "copy"
  
  input:
    tuple val(meta), path(reads), val(is_paired_end), path(fasta_amb), path(fasta_ann), path(fasta_bwt), path(fasta_pac), path(fasta_sa)
    path(ref_fa)

  output:
    tuple val(meta), path("*.sorted.bam"), path("*.bai"), val(is_paired_end), emit: regular_output
    path("${meta.id}.trimmed_reads.txt"), emit: trimmed_reads, optional: true
    
  script:
    trim_bam = "${meta.id}.trimmed"
    bed = "${params.primersBED}"

    """
    # Link reference files
    ln -s ${ref_fa} ./${fasta_amb.getSimpleName()}.fasta

    if [[ ${is_paired_end} == true ]]; then
        bwa mem ./${ref_fa} ${reads[0]} ${reads[1]} \
                -o ${meta.id}.bam -t ${params.bwa_threads} 
    else
        bwa mem ./${ref_fa} ${reads[0]} \
                -o ${meta.id}.bam -t ${params.bwa_threads} 
    fi

    # Sort and index
    samtools sort -o ${meta.id}.sorted.bam ${meta.id}.bam
    samtools index ${meta.id}.sorted.bam

    # Trim primers if bed file is provided
    if [[ "${params.primersBED}" != "null" ]]; then
        samtools ampliconclip --both-ends --hard-clip \
          --filter-len ${params.minLen} \
          -b ${bed} ${meta.id}.sorted.bam \
          -f ${meta.id}.trimmed_reads.txt > ${trim_bam}

        samtools sort -o ${trim_bam}.sorted.bam ${trim_bam}
        samtools index ${trim_bam}.sorted.bam

        mv ${meta.id}.sorted.bam ${meta.id}.raw.sorted.bam
        mv ${meta.id}.sorted.bam.bai ${meta.id}.raw.sorted.bam.bai
        mv ${trim_bam}.sorted.bam ${meta.id}.sorted.bam
        mv ${trim_bam}.sorted.bam.bai ${meta.id}.sorted.bam.bai
    fi
    """
}

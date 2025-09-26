
process align2ref{
  publishDir "${params.outDir}/${sample_id}_results/", mode: "copy"
  
  input:
    tuple val(sample_id), path(reads), val(is_paired_end), path(fasta_amb), path(fasta_ann), path(fasta_bwt), path(fasta_pac), path(fasta_sa)
    path(ref_fa)

  output:
    tuple val(sample_id), path("*.sorted.bam"), path("*.bai"), val(is_paired_end), emit: regular_output
    path("${sample_id}.trimmed_reads.txt"), emit: trimmed_reads, optional: true
    
  script:
    trim_bam = "${sample_id}.trimmed"
    bed = "${params.primersBED}"

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

    # Trim primers if bed file is provided
    if [[ "${params.primersBED}" != "null" ]]; then
        samtools ampliconclip --both-ends --hard-clip \
          --filter-len ${params.minLen} \
          -b ${bed} ${sample_id}.sorted.bam \
          -f ${sample_id}.trimmed_reads.txt > ${trim_bam}

        samtools sort -o ${trim_bam}.sorted.bam ${trim_bam}
        samtools index ${trim_bam}.sorted.bam

        mv ${sample_id}.sorted.bam ${sample_id}.raw.sorted.bam
        mv ${sample_id}.sorted.bam.bai ${sample_id}.raw.sorted.bam.bai
        mv ${trim_bam}.sorted.bam ${sample_id}.sorted.bam
        mv ${trim_bam}.sorted.bam.bai ${sample_id}.sorted.bam.bai
    fi
    """
}

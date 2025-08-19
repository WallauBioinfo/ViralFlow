process run_clair3{
    // Define the process parameters
    publishDir "${params.outDir}/${meta.id}/", mode: 'copy', overwrite: true
    tag "${meta.id}"
    container "docker://hkubal/clair3:v1.1.0"
    
    input:
        tuple val(meta), path(bam)
        path(ref)
        val(chunk_size) // 10000
        val(qual) // 30
        val(model)

    output:
        tuple val(meta), path("${meta.id}.merge_output.vcf.gz")
  
    script:
    """
    samtools faidx ${ref}
    
    run_clair3.sh \
        --enable_long_indel \
        --chunk_size=${chunk_size} \
        --haploid_sensitive \
        --no_phasing_for_fa \
        --bam_fn=${bam} \
        --ref_fn=${ref} \
        --output=./ \
        --threads=${task.cpus} \
        --platform='ont' \
        --model_path=/opt/models/${model} \
        --include_all_ctgs \
        --min_mq=${qual}
    
    mv ./merge_output.vcf.gz ${meta.id}.merge_output.vcf.gz
    """
}

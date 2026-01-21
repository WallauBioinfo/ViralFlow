
process run_bcftools {
    label "NP_basecontainer"
    // Define the process parameters
    publishDir "${params.outDir}/${meta.id}_results/", mode: 'copy', overwrite: true
    tag "${meta.id}"

    input:
        tuple val(meta), path(vcf)
        path(ref)
        val(af_threshold) // 0.51

    output:
        tuple val(meta), path("${meta.id}.filtered.vcf.gz*")
  
    script:
    """
    # normalize indels
    bcftools norm -m - -f ${ref} ${vcf} > norm.vcf
    
    # filter variants based on allele frequency
    bcftools filter -i "FORMAT/AF >= ${af_threshold}" norm.vcf -o ${meta.id}.filtered.vcf -O z
    
    # index the filtered VCF
    bgzip ${meta.id}.filtered.vcf
    tabix ${meta.id}.filtered.vcf.gz
    """
}

process run_bcftools_consensus {
    label "NP_basecontainer"
    // Define the process parameters
    publishDir "${params.outDir}/${meta.id}_results/", mode: 'copy', overwrite: true
    tag "${meta.id}"

    input:
        tuple val(meta), path(vcf), path(sorted_bam_files)
        path(ref)
        val(min_depth) // 3
    
    output:
        tuple val(meta), path("${meta.id}.consensus.fa"), path("${meta.id}.low_cov.bed"), path("${meta.id}.cov.bed")
  
    shell:
    vcf_file = "${vcf[0]}"
    sorted_bam = "${sorted_bam_files[0]}"
    '''
    # create a bed file with low coverage regions
    # this is used to mask low coverage regions in the consensus sequence
    samtools depth -J -a !{sorted_bam} > !{meta.id}.cov.bed
    awk '$3 <= int(!{min_depth}) {print $1 "\t" $2-1 "\t" $2}' !{meta.id}.cov.bed > !{meta.id}.low_cov.bed
    
    # call consensus sequence and rename it
    bcftools consensus -f !{ref} --mask !{meta.id}.low_cov.bed !{vcf_file} > !{meta.id}.consensus.fa
    sed -i -e 's/>.*/>!{meta.id}/' !{meta.id}.consensus.fa
    '''
}

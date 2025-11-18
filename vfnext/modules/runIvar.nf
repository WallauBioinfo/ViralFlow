process runIvar{
  publishDir "${params.outDir}/${sample_id}_results/", mode: "copy", pattern: "*.{fa,tsv,gz,tbi}"
  input:
    tuple val(sample_id), path(bams), val(is_paired_end)
    path(ref_fa)
  output:
    tuple val(sample_id), path("*.depth*.fa"), path("*.txt"), path("${sample_id}.tsv"), path("${sample_id}.ivar.vcf.gz"), path("${sample_id}.ivar.vcf.gz.tbi")

  script:
    sorted_bam = "${bams[0].getSimpleName()}.sorted.bam"
    d = "${params.depth}"
    """
    # IVAR STEP 1 ----------------------------------------------------------------
    samtools mpileup -aa -d 50000 --reference ${ref_fa} -a -B ${sorted_bam} | \
       ivar variants -G -p ${sample_id} -q ${params.mapping_quality} -t 0.05
    python $projectDir/bin/tsv_to_vcf.py ${sample_id}.tsv ${sample_id}.ivar.vcf ${sample_id}
    bgzip ${sample_id}.ivar.vcf
    tabix ${sample_id}.ivar.vcf.gz

    # IVAR STEP 2 ----------------------------------------------------------------
    samtools mpileup -aa -d 50000 --reference ${ref_fa} -a -B ${sorted_bam} | \
       ivar consensus -p ${sample_id} -q ${params.mapping_quality} -t 0 -m ${d} -n N -c 0.51

    # IVAR STEP 3 ----------------------------------------------------------------
    samtools mpileup -aa -d 50000 --reference ${ref_fa} -a -B ${sorted_bam} | \
       ivar consensus -p ${sample_id}.ivar060 -q ${params.mapping_quality} -t 0.60 -n N -m ${params.depth} -c 0.51
    # EDIT FILE NAMES
    mv ${sample_id}.fa ${sample_id}.depth${d}.fa
    mv ${sample_id}.ivar060.fa ${sample_id}.depth${d}.amb.fa
    sed -i -e 's/>.*/>${sample_id}/g' ${sample_id}.depth${d}.fa
    sed -i -e 's/>.*/>${sample_id}/g' ${sample_id}.depth${d}.amb.fa
    """
}

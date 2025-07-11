process runIvar{
  publishDir "${params.outDir}/${sample_id}_results/", mode: "copy", pattern: "*.{fa,tsv}"
  input:
    tuple val(sample_id), path(bams), val(is_paired_end)
    path(ref_fa)
  output:
    tuple val(sample_id), path("*.depth*.fa"), path("*.txt"), path("${sample_id}.tsv")

  script:
    sorted_bam = "${bams[0].getSimpleName()}.sorted.bam"
    d = "${params.depth}"
    """
    # IVAR STEP 1 ----------------------------------------------------------------
    samtools mpileup -aa -d 50000 --reference ${ref_fa} -a -B ${sorted_bam} | \
       ivar variants -p ${sample_id} -q ${params.mapping_quality} -t 0.05

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

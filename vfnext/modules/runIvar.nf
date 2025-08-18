process runIvar{
  tag "${meta.id}"
  publishDir "${params.outDir}/${meta.id}_results/", mode: "copy", pattern: "*.{fa,tsv}"
  input:
    tuple val(meta), path(bams), val(is_paired_end)
    path(ref_fa)
  output:
    tuple val(meta), path("*.depth*.fa"), path("*.txt"), path("${meta.id}.tsv")

  script:
    sorted_bam = "${bams[0].getSimpleName()}.sorted.bam"
    d = "${params.depth}"
    """
    # IVAR STEP 1 ----------------------------------------------------------------
    samtools mpileup -aa -d 50000 --reference ${ref_fa} -a -B ${sorted_bam} | \
       ivar variants -p ${meta.id} -q ${params.mapping_quality} -t 0.05

    # IVAR STEP 2 ----------------------------------------------------------------
    samtools mpileup -aa -d 50000 --reference ${ref_fa} -a -B ${sorted_bam} | \
       ivar consensus -p ${meta.id} -q ${params.mapping_quality} -t 0 -m ${d} -n N -c 0.51

    # IVAR STEP 3 ----------------------------------------------------------------
    samtools mpileup -aa -d 50000 --reference ${ref_fa} -a -B ${sorted_bam} | \
       ivar consensus -p ${meta.id}.ivar060 -q ${params.mapping_quality} -t 0.60 -n N -m ${params.depth} -c 0.51
    # EDIT FILE NAMES
    mv ${meta.id}.fa ${meta.id}.depth${d}.fa
    mv ${meta.id}.ivar060.fa ${meta.id}.depth${d}.amb.fa
    sed -i -e 's/>.*/>${meta.id}/g' ${meta.id}.depth${d}.fa
    sed -i -e 's/>.*/>${meta.id}/g' ${meta.id}.depth${d}.amb.fa
    """
}

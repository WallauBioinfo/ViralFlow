process genFaIdx {
    /*
    * Indexes reference fasta file using bwa.
    */
    //publishDir "${workflow.outputDir}/"
    label "singlethread"
    input:
        path(reference_fasta)

    output:
        path("${reference_fasta}*")

    script:
        """
        samtools faidx ${reference_fasta}
        """
}

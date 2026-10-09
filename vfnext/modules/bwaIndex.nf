params.bwa='bwa'

process indexReferenceBWA {
    /**
    * Indexes reference fasta file using bwa.
    */
    //publishDir "${params.outDir}/"
    label "singlethread"

    input:
        path(refFa)

    output:
        path("${refFa}*")

    script:
        bwa=params.bwa

        """
        ${bwa} index -a bwtsw -p ${refFa} ${refFa}
        """
}

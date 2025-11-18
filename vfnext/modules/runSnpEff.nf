process runSnpEff{
  errorStrategy 'ignore'
  publishDir "${params.outDir}/${sample_id}_results/", mode: "copy"

        input:
        val(genome_code)
        path(entry_found_file) // here to assure runSnpEff will not start before DB was checked
        tuple val(sample_id), path(consensus), path(ivar_txt), path(mut_tsv),  path(vcf_file), path(vcf_index)
        output:
        tuple val(sample_id), path("*.vcf"), path("${sample_id}_snpEff_summary.html"), path("snpEff_genes.txt")

        script:
        if (params.virus == "custom")
        """        
        snpEff ann -Xmx4g \
                -v ${genome_code} ${vcf_file} > ${sample_id}.ann.vcf 2> snpEff.log
        # add sample id to htmls
        mv snpEff_summary.html ${sample_id}_snpEff_summary.html
        """
        else
        """
        snpEff ann -Xmx4g \
                ${genome_code} ${vcf_file} > ${sample_id}.ann.vcf 2> snpEff.log
        # add sample id to htmls
        mv snpEff_summary.html ${sample_id}_snpEff_summary.html
        """
}

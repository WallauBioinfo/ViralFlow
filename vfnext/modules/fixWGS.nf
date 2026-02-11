process fixWGS {
  tag "${meta.id}"
  label "singlethread"
  errorStrategy 'ignore'
  publishDir "${params.outDir}/${meta.id}_results/", mode : "copy"

  input:
     tuple val(sample_id), path(wgs), path(metrics), path(consensus), path(ivar_txt), path(mut_tsv), path(vcf_file), path(vcf_index)
     //temporary solution, no need for ivar_txt and mut_tsv vcf_file and vcf_index

  output:
     path("${meta.id}.metrics.genome.tsv")

  script:
     consensus_fa = "${meta.id}.depth${params.depth}.fa"
     """
     #!/usr/bin/env python
     # ----- import libraries -------------------------------------------------
     import pandas as pd
     from Bio import SeqIO
     import sys
     # ----- Function ---------------------------------------------------------

     def computeCoverage(consensus):
        '''
        (N - total bases) / total_bases
        '''
        for record in SeqIO.parse(consensus, "fasta"):
           seq = record.seq
        total_N = sum([1 for i in seq if i == "N"])
        total_bases = len(seq)
        try:
            assert(total_bases > 0)      
        except(AssertionError):
            print("WARN: No sequence at ${consensus_fa}")
            return 0
        
        return (total_bases - total_N) / total_bases
     # ------------------------------------------------------------------------
     # compute coverage
     exact_coverage = computeCoverage("${consensus_fa}")
     # add to wgs
     wgs_tsv = pd.read_csv(f"${wgs}", comment="#", sep='\t', skiprows=2)
     wgs_tsv = wgs_tsv[:1]
     wgs_tsv['EXACT_COV'] = exact_coverage
     # write to a new tsv file
     with open(f"${meta.id}.metrics.genome.tsv",'w') as wgs_out:
         wgs_tsv.to_csv(wgs_out, sep ='\t', index = False)
     """
}

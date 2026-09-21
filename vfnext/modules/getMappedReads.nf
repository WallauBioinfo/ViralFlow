process getMappedReads{
  tag "${meta.id}"
  publishDir { "${params.outDir}/${meta.id}_results/" }, mode: "copy"
  label "singlethread"
  input:
    tuple val(meta), path(bam), val(is_paired_end)

  output:
    tuple val(meta), path("*.mapped.*.fq.gz")
  script:
    """
    if [[ ${is_paired_end}  == true ]]; then
      samtools sort -n ${bam} | \
      samtools fastq -F 4 -1 ${meta.id}.mapped.R1.fq.gz -2 ${meta.id}.mapped.R2.fq.gz
    else
      samtools sort -n ${bam} | \
      samtools fastq -F 4 -0 ${meta.id}.mapped.SE.fq.gz
    fi
    """
}

/*
// --- DOCUMENTATION ----
This process was designed to get fastqs containing only the mapped reads

1 - "samtools sort -n" organize reads by name (the "sorted.bam" is organized by quality)
2 - "samtools fastq" filter the mapped reads and write it as fastq

On the single-end branch the reads are collected with "-0", not "-s". The two
are not interchangeable:

  -s  writes SINGLETONS - reads that carry the paired flag but whose mate is
      gone. Single-end reads are not flagged as paired, so they are never
      singletons and "-s" receives nothing; samtools sends them to stdout
      instead. Asking for "-s out.fq.gz" on single-end data therefore publishes
      an empty archive and drops every read into .command.out. Verified with
      samtools 1.21: 3 mapped single-end reads in, a 28-byte empty gzip out.
  -0  writes reads with both or neither of the READ1/READ2 flags set, which is
      exactly what a single-end read is.

"-0" also compresses according to the file extension, so no separate gzip step
is needed. Shell redirection does not: "> out.fq.gz" writes plain text under a
.gz name, which is a trap rather than a shortcut.

Known limitation, inherited rather than introduced: the paired branch has no
"-s", so a read whose mate was removed by "-F 4" is written to the R1 file
instead of being set aside. R1 and R2 then hold different numbers of records
and no longer line up. Fixing that means adding "-s" here, which adds an output
file, so it is left for the ILLUMINA work - see TODO.md.
*/

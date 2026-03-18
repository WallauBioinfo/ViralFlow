#!/bin/bash -ue
# Link reference files
if [[ ! -f ./NC_045512.fasta ]]; then
    ln -s NC_045512.2.fa ./NC_045512.fasta
fi

if [[ true == true ]]; then
    bwa mem ./NC_045512.2.fa Cneg.R1.fq.gz Cneg.R2.fq.gz                 -o Cneg.bam -t 1 
else
    bwa mem ./NC_045512.2.fa Cneg.R1.fq.gz                 -o Cneg.bam -t 1 
fi

# Sort and index
samtools sort -o Cneg.sorted.bam Cneg.bam
samtools index Cneg.sorted.bam

# Trim primers if bed file is provided
if [[ "null" != "null" ]]; then
    samtools ampliconclip --both-ends --hard-clip           --filter-len 75           -b null Cneg.sorted.bam           -f Cneg.trimmed_reads.txt > Cneg.trimmed

    samtools sort -o Cneg.trimmed.sorted.bam Cneg.trimmed
    samtools index Cneg.trimmed.sorted.bam

    mv Cneg.sorted.bam Cneg.raw.sorted.bam
    mv Cneg.sorted.bam.bai Cneg.raw.sorted.bam.bai
    mv Cneg.trimmed.sorted.bam Cneg.sorted.bam
    mv Cneg.trimmed.sorted.bam.bai Cneg.sorted.bam.bai
fi

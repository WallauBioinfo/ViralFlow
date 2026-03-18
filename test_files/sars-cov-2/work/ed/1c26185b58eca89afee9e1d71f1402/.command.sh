#!/bin/bash -ue
# Link reference files
if [[ ! -f ./NC_045512.fasta ]]; then
    ln -s NC_045512.2.fa ./NC_045512.fasta
fi

if [[ false == true ]]; then
    bwa mem ./NC_045512.2.fa ART3.SE.fq.gz null                 -o ART3.bam -t 1 
else
    bwa mem ./NC_045512.2.fa ART3.SE.fq.gz                 -o ART3.bam -t 1 
fi

# Sort and index
samtools sort -o ART3.sorted.bam ART3.bam
samtools index ART3.sorted.bam

# Trim primers if bed file is provided
if [[ "null" != "null" ]]; then
    samtools ampliconclip --both-ends --hard-clip           --filter-len 75           -b null ART3.sorted.bam           -f ART3.trimmed_reads.txt > ART3.trimmed

    samtools sort -o ART3.trimmed.sorted.bam ART3.trimmed
    samtools index ART3.trimmed.sorted.bam

    mv ART3.sorted.bam ART3.raw.sorted.bam
    mv ART3.sorted.bam.bai ART3.raw.sorted.bam.bai
    mv ART3.trimmed.sorted.bam ART3.sorted.bam
    mv ART3.trimmed.sorted.bam.bai ART3.sorted.bam.bai
fi

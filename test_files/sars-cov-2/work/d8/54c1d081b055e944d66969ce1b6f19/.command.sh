#!/bin/bash -ue
# Link reference files
if [[ ! -f ./NC_045512.fasta ]]; then
    ln -s NC_045512.2.fa ./NC_045512.fasta
fi

if [[ true == true ]]; then
    bwa mem ./NC_045512.2.fa ART1.R1.fq.gz ART1.R2.fq.gz                 -o ART1.bam -t 1 
else
    bwa mem ./NC_045512.2.fa ART1.R1.fq.gz                 -o ART1.bam -t 1 
fi

# Sort and index
samtools sort -o ART1.sorted.bam ART1.bam
samtools index ART1.sorted.bam

# Trim primers if bed file is provided
if [[ "null" != "null" ]]; then
    samtools ampliconclip --both-ends --hard-clip           --filter-len 75           -b null ART1.sorted.bam           -f ART1.trimmed_reads.txt > ART1.trimmed

    samtools sort -o ART1.trimmed.sorted.bam ART1.trimmed
    samtools index ART1.trimmed.sorted.bam

    mv ART1.sorted.bam ART1.raw.sorted.bam
    mv ART1.sorted.bam.bai ART1.raw.sorted.bam.bai
    mv ART1.trimmed.sorted.bam ART1.sorted.bam
    mv ART1.trimmed.sorted.bam.bai ART1.sorted.bam.bai
fi

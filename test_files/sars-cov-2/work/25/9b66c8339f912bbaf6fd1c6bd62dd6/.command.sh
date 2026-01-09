#!/bin/bash -ue
# Link reference files
if [[ ! -f ./NC_045512.fasta ]]; then
    ln -s NC_045512.2.fa ./NC_045512.fasta
fi

if [[ true == true ]]; then
    bwa mem ./NC_045512.2.fa ART2.R1.fq.gz ART2.R2.fq.gz                 -o ART2.bam -t 1 
else
    bwa mem ./NC_045512.2.fa ART2.R1.fq.gz                 -o ART2.bam -t 1 
fi

# Sort and index
samtools sort -o ART2.sorted.bam ART2.bam
samtools index ART2.sorted.bam

# Trim primers if bed file is provided
if [[ "null" != "null" ]]; then
    samtools ampliconclip --both-ends --hard-clip           --filter-len 75           -b null ART2.sorted.bam           -f ART2.trimmed_reads.txt > ART2.trimmed

    samtools sort -o ART2.trimmed.sorted.bam ART2.trimmed
    samtools index ART2.trimmed.sorted.bam

    mv ART2.sorted.bam ART2.raw.sorted.bam
    mv ART2.sorted.bam.bai ART2.raw.sorted.bam.bai
    mv ART2.trimmed.sorted.bam ART2.sorted.bam
    mv ART2.trimmed.sorted.bam.bai ART2.sorted.bam.bai
fi

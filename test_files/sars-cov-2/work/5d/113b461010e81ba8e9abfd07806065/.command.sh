#!/bin/bash -ue
samtools view -b -f 4 ART2.sorted.bam > unmapped.bam
if [[ true  == true ]]; then
  samtools sort -n unmapped.bam |       samtools fastq -f 4 -1 ART2.unmapped.R1.fq.gz -2 ART2.unmapped.R2.fq.gz
else
  samtools sort -n unmapped.bam |       samtools fastq -f 4 -s ART2.unmapped.SE.fq.gz
fi

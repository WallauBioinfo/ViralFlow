#!/bin/bash -ue
samtools view -b -f 4 ART1.sorted.bam > unmapped.bam
if [[ true  == true ]]; then
  samtools sort -n unmapped.bam |       samtools fastq -f 4 -1 ART1.unmapped.R1.fq.gz -2 ART1.unmapped.R2.fq.gz
else
  samtools sort -n unmapped.bam |       samtools fastq -f 4 -s ART1.unmapped.SE.fq.gz
fi

#!/bin/bash -ue
samtools view -b -f 4 ART3.sorted.bam > unmapped.bam
if [[ false  == true ]]; then
  samtools sort -n unmapped.bam |       samtools fastq -f 4 -1 ART3.unmapped.R1.fq.gz -2 ART3.unmapped.R2.fq.gz
else
  samtools sort -n unmapped.bam |       samtools fastq -f 4 -s ART3.unmapped.SE.fq.gz
fi

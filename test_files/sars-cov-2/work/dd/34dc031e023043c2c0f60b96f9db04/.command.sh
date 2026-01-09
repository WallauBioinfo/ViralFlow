#!/bin/bash -ue
if [[ true  == true ]]; then
  samtools sort -n ART1.sorted.bam |       samtools fastq -F 4 -1 ART1.mapped.R1.fq.gz -2 ART1.mapped.R2.fq.gz
else
  samtools sort -n ART1.sorted.bam |       samtools fastq -F 4 -s ART1.mapped.SE.fq.gz 
fi

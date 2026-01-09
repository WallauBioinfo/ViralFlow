#!/bin/bash -ue
if [[ false  == true ]]; then
  samtools sort -n ART3.sorted.bam |       samtools fastq -F 4 -1 ART3.mapped.R1.fq.gz -2 ART3.mapped.R2.fq.gz
else
  samtools sort -n ART3.sorted.bam |       samtools fastq -F 4 -s ART3.mapped.SE.fq.gz 
fi

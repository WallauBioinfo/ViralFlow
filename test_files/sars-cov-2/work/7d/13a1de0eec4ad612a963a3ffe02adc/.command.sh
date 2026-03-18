#!/bin/bash -ue
if [[ true  == true ]]; then
  samtools sort -n ART2.sorted.bam |       samtools fastq -F 4 -1 ART2.mapped.R1.fq.gz -2 ART2.mapped.R2.fq.gz
else
  samtools sort -n ART2.sorted.bam |       samtools fastq -F 4 -s ART2.mapped.SE.fq.gz 
fi

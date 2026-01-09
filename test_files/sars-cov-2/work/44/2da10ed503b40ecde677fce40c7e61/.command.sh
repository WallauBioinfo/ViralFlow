#!/bin/bash -ue
if [[ true  == true ]]; then
  samtools sort -n Cneg.sorted.bam |       samtools fastq -F 4 -1 Cneg.mapped.R1.fq.gz -2 Cneg.mapped.R2.fq.gz
else
  samtools sort -n Cneg.sorted.bam |       samtools fastq -F 4 -s Cneg.mapped.SE.fq.gz 
fi

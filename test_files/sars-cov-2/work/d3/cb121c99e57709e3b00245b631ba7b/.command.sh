#!/bin/bash -ue
samtools view -b -f 4 Cneg.sorted.bam > unmapped.bam
if [[ true  == true ]]; then
  samtools sort -n unmapped.bam |       samtools fastq -f 4 -1 Cneg.unmapped.R1.fq.gz -2 Cneg.unmapped.R2.fq.gz
else
  samtools sort -n unmapped.bam |       samtools fastq -f 4 -s Cneg.unmapped.SE.fq.gz
fi

#!/bin/bash -ue
freebayes -p 1 --reference-quality 30,30             -f NC_045512.2.fa ART1.sorted.bam > ART1.vcf
snpEff ann -Xmx4g             NC_045512.2 ART1.vcf > ART1.ann.vcf
# add sample id to htmls
mv snpEff_summary.html ART1_snpEff_summary.html

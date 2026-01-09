#!/bin/bash -ue
freebayes -p 1 --reference-quality 30,30             -f NC_045512.2.fa ART2.sorted.bam > ART2.vcf
snpEff ann -Xmx4g             NC_045512.2 ART2.vcf > ART2.ann.vcf
# add sample id to htmls
mv snpEff_summary.html ART2_snpEff_summary.html

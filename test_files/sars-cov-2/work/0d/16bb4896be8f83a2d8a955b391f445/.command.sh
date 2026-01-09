#!/bin/bash -ue
freebayes -p 1 --reference-quality 30,30             -f NC_045512.2.fa ART3.sorted.bam > ART3.vcf
snpEff ann -Xmx4g             NC_045512.2 ART3.vcf > ART3.ann.vcf
# add sample id to htmls
mv snpEff_summary.html ART3_snpEff_summary.html

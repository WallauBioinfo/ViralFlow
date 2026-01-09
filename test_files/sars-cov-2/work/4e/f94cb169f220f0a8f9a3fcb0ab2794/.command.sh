#!/bin/bash -ue
freebayes -p 1 --reference-quality 30,30             -f NC_045512.2.fa Cneg.sorted.bam > Cneg.vcf
snpEff ann -Xmx4g             NC_045512.2 Cneg.vcf > Cneg.ann.vcf
# add sample id to htmls
mv snpEff_summary.html Cneg_snpEff_summary.html

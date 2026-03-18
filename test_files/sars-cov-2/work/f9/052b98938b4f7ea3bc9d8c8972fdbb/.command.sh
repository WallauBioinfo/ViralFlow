#!/usr/bin/env python
# ----- import libraries -------------------------------------------------
import pysam
import subprocess

#Loading the BAM file
bam_file = pysam.AlignmentFile('ART1.sorted.bam', 'rb')

#Checking the number of mapped reads
n_mapped_reads = bam_file.count()

if n_mapped_reads > 0:

  references = bam_file.references
  if len(references) == 0:
    raise ValueError("No references found in BAM header.")

  genomecode = references[0]

  subprocess.run([f"bamdash -b ART1.sorted.bam -r {genomecode} -c 25 -e svg"], shell=True)
  subprocess.run([f"mv {genomecode}_plot.svg ART1_coveragePlot.svg"], shell=True)

  subprocess.run([f"bamdash -b ART1.sorted.bam -r {genomecode} -c 25 -e png"], shell=True)
  subprocess.run([f"mv {genomecode}_plot.png ART1_coveragePlot.png"], shell=True)

  subprocess.run([f"bamdash -b ART1.sorted.bam -r {genomecode} -c 25"], shell=True)
  subprocess.run([f"mv {genomecode}_plot.html ART1_coveragePlot.html"], shell=True)

else:
  result = "No mapped reads were found in the sorted BAM file for sample ART1. The coverage plot will not be generated for it."
  with open('coveragePlot_result.txt', 'w') as f:
    f.write(result)

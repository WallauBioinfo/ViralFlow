#!/bin/bash -ue
# IVAR STEP 1 ----------------------------------------------------------------
samtools mpileup -aa -d 50000 --reference NC_045512.2.fa -a -B Cneg.sorted.bam |        ivar variants -p Cneg -q 30 -t 0.05

# IVAR STEP 2 ----------------------------------------------------------------
samtools mpileup -aa -d 50000 --reference NC_045512.2.fa -a -B Cneg.sorted.bam |        ivar consensus -p Cneg -q 30 -t 0 -m 25 -n N -c 0.51

# IVAR STEP 3 ----------------------------------------------------------------
samtools mpileup -aa -d 50000 --reference NC_045512.2.fa -a -B Cneg.sorted.bam |        ivar consensus -p Cneg.ivar060 -q 30 -t 0.60 -n N -m 25 -c 0.51
# EDIT FILE NAMES
mv Cneg.fa Cneg.depth25.fa
mv Cneg.ivar060.fa Cneg.depth25.amb.fa
sed -i -e 's/>.*/>Cneg/g' Cneg.depth25.fa
sed -i -e 's/>.*/>Cneg/g' Cneg.depth25.amb.fa
